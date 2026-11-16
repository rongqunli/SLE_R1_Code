if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(c("data.table", "limma", "matrixStats"),
                 "longitudinal gene associations")

script_path <- script_file()
expression_file <- file.path(cfg$paths$processed, "01_gene_expression.rds")
clinical_file <- file.path(
  cfg$paths$results, "03_longitudinal", "longitudinal_analysis_data.tsv"
)
output_dir <- file.path(cfg$paths$results, "05_gene_associations")
object_file <- file.path(cfg$paths$processed, "03_gene_associations.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- c(
  object_file,
  file.path(output_dir, c(
    "longitudinal_gene_results.tsv.gz", "leading_gene_results.tsv",
    "gene_model_summary.tsv", "phenotype_centering_audit.tsv",
    "correlation_estimation_audit.tsv"
  ))
)
signature <- build_step_signature(
  input_files = c(expression_file, clinical_file, cfg$files$programs),
  code_files = c(script_path, file.path(PIPELINE_ROOT, "scripts", "_common.R"),
                 file.path(PIPELINE_ROOT, "config", "config.R")),
  extra = list(method = "limma_duplicateCorrelation_v1", n_model_families = 12L)
)
if (!has_flag("--force") && step_is_complete(output_dir, signature, outputs)) {
  log_message("Step 05 gene associations are already complete.")
  quit(save = "no", status = 0L)
}

safe_write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid(), ".gz")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  data.table::fwrite(
    x, temporary, sep = "\t", quote = FALSE, na = "NA", compress = "gzip"
  )
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(temporary, path)) stop("Could not write compressed TSV: ", path)
  invisible(path)
}

parse_numeric_field <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "not_available", "not_applicable", "Data Not Available", "unknown")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

center_patient_value <- function(data, value_column, prefix) {
  between_column <- paste0(prefix, "_between")
  within_column <- paste0(prefix, "_within")
  data[, (between_column) := {
    observed <- get(value_column)
    if (all(is.na(observed))) NA_real_ else mean(observed, na.rm = TRUE)
  }, by = patient_id]
  data[, (within_column) := get(value_column) - get(between_column)]
  invisible(data)
}

fit_gene_family <- function(
    expression, data, axis, family, rhs, target_term, predictor,
    correlation, formula_label = NULL) {
  design_formula <- stats::as.formula(paste("~", rhs))
  variables <- all.vars(design_formula)
  complete <- stats::complete.cases(data[, ..variables]) & !is.na(data$patient_id)
  model_data <- droplevels(as.data.frame(data[complete]))
  model_expression <- expression[, complete, drop = FALSE]
  design <- stats::model.matrix(design_formula, data = model_data)
  if (qr(design)$rank != ncol(design)) {
    stop("Rank-deficient design for model family: ", family)
  }
  if (!target_term %in% colnames(design)) {
    stop("Target term missing from design for ", family, ": ", target_term)
  }
  target_index <- match(target_term, colnames(design))
  block <- factor(model_data$patient_id)
  fit <- limma::lmFit(
    model_expression, design = design, block = block, correlation = correlation
  )
  fit <- limma::eBayes(fit, trend = TRUE)

  estimate <- unname(fit$coefficients[, target_index])
  standard_error <- sqrt(fit$s2.post) * fit$stdev.unscaled[, target_index]
  t_statistic <- unname(fit$t[, target_index])
  p_value <- unname(fit$p.value[, target_index])
  degrees_freedom <- fit$df.total
  critical_value <- stats::qt(0.975, df = degrees_freedom)
  predictor_values <- predictor[complete]
  predictor_sd <- stats::sd(as.numeric(predictor_values), na.rm = TRUE)
  if (!is.finite(predictor_sd) || predictor_sd <= 0) {
    stop("Predictor has no finite variation for model family: ", family)
  }
  formula_text <- formula_label %||% paste0(
    "expression ~ ", rhs, " + patient block (limma duplicateCorrelation)"
  )

  data.table::data.table(
    axis = axis,
    model_family = family,
    gene = rownames(model_expression),
    term = target_term,
    estimate = estimate,
    std_error = standard_error,
    conf_low = estimate - critical_value * standard_error,
    conf_high = estimate + critical_value * standard_error,
    t_statistic = t_statistic,
    p_value = p_value,
    fdr = stats::p.adjust(p_value, method = "BH"),
    average_expression = rowMeans(model_expression),
    effect_per_sd = estimate * predictor_sd,
    predictor_sd = predictor_sd,
    n_patients = data.table::uniqueN(model_data$patient_id),
    n_visits = nrow(model_data),
    consensus_correlation = correlation,
    formula = formula_text
  )
}

select_correlation_genes <- function(expression, n_genes = 500L) {
  variances <- matrixStats::rowVars(expression)
  ordered <- order(variances, na.last = NA)
  positions <- unique(round(seq(1, length(ordered), length.out = n_genes)))
  ordered[positions]
}

estimate_consensus_correlation <- function(
    expression, data, rhs, sample_set, n_genes = 500L) {
  design_formula <- stats::as.formula(paste("~", rhs))
  variables <- all.vars(design_formula)
  complete <- stats::complete.cases(data[, ..variables]) & !is.na(data$patient_id)
  model_data <- droplevels(as.data.frame(data[complete]))
  model_expression <- expression[, complete, drop = FALSE]
  design <- stats::model.matrix(design_formula, data = model_data)
  if (qr(design)$rank != ncol(design)) {
    stop("Rank-deficient correlation design for sample set: ", sample_set)
  }
  selected <- select_correlation_genes(model_expression, n_genes)
  elapsed <- system.time({
    correlation <- limma::duplicateCorrelation(
      model_expression[selected, , drop = FALSE], design = design,
      block = factor(model_data$patient_id)
    )$consensus.correlation
  })[["elapsed"]]
  if (!is.finite(correlation)) correlation <- 0
  list(
    value = correlation,
    audit = data.table::data.table(
      sample_set = sample_set,
      n_samples = nrow(model_data),
      n_patients = data.table::uniqueN(model_data$patient_id),
      n_genes_used = length(selected),
      consensus_correlation = correlation,
      elapsed_seconds = as.numeric(elapsed),
      formula = paste0("correlation ~ ", rhs)
    )
  )
}

expression_object <- readRDS(expression_file)
expression <- expression_object$expression
if (!is.matrix(expression) || nrow(expression) != 19307L || ncol(expression) != 972L) {
  stop("Unexpected gene-expression object dimensions.")
}
metadata <- data.table::as.data.table(expression_object$metadata)
clinical <- data.table::fread(clinical_file, sep = "\t", na.strings = "NA")
if (anyDuplicated(metadata$sample_id) || anyDuplicated(clinical$sample_id)) {
  stop("Gene association inputs contain duplicate sample IDs.")
}

sle_ids <- clinical$sample_id
sle_columns <- match(sle_ids, colnames(expression))
if (anyNA(sle_columns)) stop("SLE clinical samples are missing from expression matrix.")
sle_expression <- expression[, sle_columns, drop = FALSE]
sle_data <- data.table::copy(clinical)
sle_data[, `:=`(
  patient_id = factor(patient_id),
  batch = factor(as.character(batch)),
  renal_numeric = as.numeric(renal_activity_numeric)
)]
center_patient_value(sle_data, "renal_numeric", "renal")

composition_covariates <- c(
  "NK", "T_CELL", "NEUTROPHIL", "MONOCYTE", "ERYTHROID", "PLATELET", "B_CELL"
)
if (!all(composition_covariates %in% names(sle_data))) {
  stop("Missing composition covariates in longitudinal analysis data.")
}
composition_rhs <- paste(composition_covariates, collapse = " + ")

all_columns <- match(colnames(expression), metadata$sample_id)
if (anyNA(all_columns)) stop("Expression metadata sample order is incomplete.")
all_data <- data.table::copy(metadata[all_columns])
all_data[, `:=`(
  patient_id = factor(patient_id),
  disease_state = factor(disease_state, levels = c("Healthy", "SLE")),
  age_numeric = parse_numeric_field(age),
  sex = factor(as.character(sex)),
  batch = factor(as.character(batch))
)]

correlation_specs <- list(
  sle_all = list(
    expression = sle_expression, data = sle_data,
    rhs = paste(
      "SLEDAI_within + SLEDAI_between + renal_within + renal_between +",
      "timepoint_scaled + batch +", composition_rhs
    )
  ),
  sle_c3 = list(
    expression = sle_expression, data = sle_data,
    rhs = paste(
      "C3_within + C3_between + SLEDAI_within + SLEDAI_between +",
      "timepoint_scaled + batch +", composition_rhs
    )
  ),
  sle_c4 = list(
    expression = sle_expression, data = sle_data,
    rhs = paste(
      "C4_within + C4_between + SLEDAI_within + SLEDAI_between +",
      "timepoint_scaled + batch +", composition_rhs
    )
  ),
  all_samples = list(
    expression = expression, data = all_data,
    rhs = "disease_state + age_numeric + sex + batch"
  )
)
correlation_cache <- list()
correlation_audit_list <- list()
for (sample_set in names(correlation_specs)) {
  definition <- correlation_specs[[sample_set]]
  log_message("Estimating patient correlation for sample set ", sample_set, ".")
  estimated <- estimate_consensus_correlation(
    definition$expression, definition$data, definition$rhs, sample_set
  )
  correlation_cache[[sample_set]] <- estimated$value
  correlation_audit_list[[sample_set]] <- estimated$audit
}
correlation_audit <- data.table::rbindlist(correlation_audit_list, use.names = TRUE)

specifications <- list(
  list(
    axis = "general_activity", family = "activity_base",
    rhs = "SLEDAI_within + SLEDAI_between + timepoint_scaled + batch",
    term = "SLEDAI_within", predictor = sle_data$SLEDAI_within,
    correlation_set = "sle_all"
  ),
  list(
    axis = "general_activity", family = "activity_composition_adjusted",
    rhs = paste(
      "SLEDAI_within + SLEDAI_between + timepoint_scaled + batch +",
      composition_rhs
    ),
    term = "SLEDAI_within", predictor = sle_data$SLEDAI_within,
    correlation_set = "sle_all"
  ),
  list(
    axis = "renal_complement", family = "renal_base",
    rhs = "renal_within + renal_between + timepoint_scaled + batch",
    term = "renal_within", predictor = sle_data$renal_within,
    correlation_set = "sle_all"
  ),
  list(
    axis = "renal_complement", family = "renal_sledai_adjusted",
    rhs = paste(
      "renal_within + renal_between + SLEDAI_within + SLEDAI_between +",
      "timepoint_scaled + batch"
    ),
    term = "renal_within", predictor = sle_data$renal_within,
    correlation_set = "sle_all"
  ),
  list(
    axis = "renal_complement", family = "renal_composition_adjusted",
    rhs = paste(
      "renal_within + renal_between + timepoint_scaled + batch +",
      composition_rhs
    ),
    term = "renal_within", predictor = sle_data$renal_within,
    correlation_set = "sle_all"
  ),
  list(
    axis = "renal_complement", family = "c3_base",
    rhs = "C3_within + C3_between + timepoint_scaled + batch",
    term = "C3_within", predictor = sle_data$C3_within,
    correlation_set = "sle_c3"
  ),
  list(
    axis = "renal_complement", family = "c3_sledai_adjusted",
    rhs = paste(
      "C3_within + C3_between + SLEDAI_within + SLEDAI_between +",
      "timepoint_scaled + batch"
    ),
    term = "C3_within", predictor = sle_data$C3_within,
    correlation_set = "sle_c3"
  ),
  list(
    axis = "renal_complement", family = "c3_composition_adjusted",
    rhs = paste(
      "C3_within + C3_between + timepoint_scaled + batch +",
      composition_rhs
    ),
    term = "C3_within", predictor = sle_data$C3_within,
    correlation_set = "sle_c3"
  ),
  list(
    axis = "renal_complement", family = "c4_base",
    rhs = "C4_within + C4_between + timepoint_scaled + batch",
    term = "C4_within", predictor = sle_data$C4_within,
    correlation_set = "sle_c4"
  ),
  list(
    axis = "renal_complement", family = "c4_sledai_adjusted",
    rhs = paste(
      "C4_within + C4_between + SLEDAI_within + SLEDAI_between +",
      "timepoint_scaled + batch"
    ),
    term = "C4_within", predictor = sle_data$C4_within,
    correlation_set = "sle_c4"
  ),
  list(
    axis = "renal_complement", family = "c4_composition_adjusted",
    rhs = paste(
      "C4_within + C4_between + timepoint_scaled + batch +",
      composition_rhs
    ),
    term = "C4_within", predictor = sle_data$C4_within,
    correlation_set = "sle_c4"
  )
)

result_list <- vector("list", length(specifications) + 1L)
for (index in seq_along(specifications)) {
  specification <- specifications[[index]]
  log_message("Fitting gene model family ", specification$family, ".")
  result_list[[index]] <- fit_gene_family(
    sle_expression, sle_data, specification$axis, specification$family,
    specification$rhs, specification$term, specification$predictor,
    correlation_cache[[specification$correlation_set]]
  )
}

log_message("Fitting auxiliary SLE versus healthy reference model.")
result_list[[length(result_list)]] <- fit_gene_family(
  expression, all_data, "healthy_reference", "healthy_reference",
  "disease_state + age_numeric + sex + batch", "disease_stateSLE",
  as.numeric(all_data$disease_state == "SLE"),
  correlation_cache[["all_samples"]]
)

results <- data.table::rbindlist(result_list, use.names = TRUE, fill = TRUE)
if (any(!is.finite(results$estimate)) || any(!is.finite(results$p_value))) {
  stop("Gene model results contain non-finite estimates or P values.")
}

model_summary <- results[, .(
  axis = unique(axis),
  term = unique(term),
  n_genes = .N,
  n_fdr_05 = sum(fdr < cfg$longitudinal$fdr),
  n_positive_fdr_05 = sum(fdr < cfg$longitudinal$fdr & estimate > 0),
  n_negative_fdr_05 = sum(fdr < cfg$longitudinal$fdr & estimate < 0),
  minimum_fdr = min(fdr),
  top_gene = gene[which.min(fdr)],
  n_patients = unique(n_patients),
  n_visits = unique(n_visits),
  consensus_correlation = unique(consensus_correlation),
  formula = unique(formula)
), by = model_family]

leading_ranked <- results[order(p_value), head(.SD, 200L), by = model_family]
leading <- unique(data.table::rbindlist(list(
  results[fdr < cfg$longitudinal$fdr], leading_ranked
), use.names = TRUE), by = c("model_family", "gene"))
data.table::setorder(leading, model_family, p_value, gene)

audit_centering <- function(data, phenotype, raw_column, within_column) {
  observed <- data[is.finite(get(within_column))]
  patient_values <- observed[, .(
    changed = diff(range(get(raw_column))) > 0,
    within_sum = sum(get(within_column))
  ), by = patient_id]
  data.table::data.table(
    phenotype = phenotype,
    n_visits = nrow(observed),
    n_patients = nrow(patient_values),
    n_changing_patients = sum(patient_values$changed),
    max_abs_patient_within_sum = max(abs(patient_values$within_sum))
  )
}
centering_audit <- data.table::rbindlist(list(
  audit_centering(sle_data, "SLEDAI", "SLEDAI_numeric", "SLEDAI_within"),
  audit_centering(sle_data, "renal", "renal_numeric", "renal_within"),
  audit_centering(sle_data, "C3", "C3_numeric", "C3_within"),
  audit_centering(sle_data, "C4", "C4_numeric", "C4_within")
), use.names = TRUE)

safe_write_tsv_gz(results, outputs[[2]])
safe_write_tsv(leading, outputs[[3]])
safe_write_tsv(model_summary, outputs[[4]])
safe_write_tsv(centering_audit, outputs[[5]])
safe_write_tsv(correlation_audit, outputs[[6]])
atomic_save_rds(list(
  results = as.data.frame(results),
  model_summary = as.data.frame(model_summary),
  centering_audit = as.data.frame(centering_audit),
  correlation_audit = as.data.frame(correlation_audit),
  analysis_data = as.data.frame(sle_data)
), object_file, compress = TRUE)
write_session_info(file.path(cfg$paths$logs, "05_gene_associations_sessionInfo.txt"))
mark_step_complete(output_dir, signature, outputs)
log_message("Gene associations complete: ", nrow(results), " gene-model results.")
