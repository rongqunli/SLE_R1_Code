if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[[1]]), winslash = "/", mustWork = TRUE)), "_common.R"))
analysis_exclusions_file <- file.path(PIPELINE_ROOT, "config", "analysis_exclusions.R")
analysis_exclusion_helpers_file <- file.path(PIPELINE_ROOT, "scripts", "_analysis_exclusion_helpers.R")
source(analysis_exclusions_file, local = FALSE)
source(analysis_exclusion_helpers_file, local = FALSE)

require_packages(
  c("SeuratObject", "Matrix", "data.table", "edgeR", "ggplot2"),
  "donor-level pseudobulk and differential abundance"
)

input_file <- file.path(cfg$paths$processed, "02_annotated_seurat.rds")
pb_dir <- file.path(cfg$paths$results, "05_pseudobulk")
da_dir <- file.path(cfg$paths$results, "05_abundance")
checkpoint_dir <- file.path(cfg$paths$processed, "checkpoints", "05_pseudobulk")
aggregation_cache_path <- file.path(checkpoint_dir, "pseudobulk_aggregations.rds")
all_genes_file <- file.path(pb_dir, "pseudobulk_all_genes.rds")
significant_file <- file.path(pb_dir, "pseudobulk_significant.tsv")
pseudobulk_metadata_file <- file.path(pb_dir, "pseudobulk_metadata.tsv")
replication_file <- file.path(pb_dir, "child_adult_direction_replication.tsv")
activity_audit_file <- file.path(pb_dir, "sledai_activity_group_audit.tsv")
cell_exclusion_file <- file.path(pb_dir, "analysis_cell_exclusions.tsv")
cell_exclusion_summary_file <- file.path(pb_dir, "analysis_cell_exclusion_summary.tsv")
abundance_results_file <- file.path(da_dir, "abundance_results.tsv")
donor_proportions_file <- file.path(da_dir, "donor_cell_proportions.tsv")
abundance_plot_file <- file.path(da_dir, "child_donor_cell_proportions.pdf")

if (!file.exists(input_file)) {
  stop("Annotated object not found: ", input_file, ". Run step 04_annotate first.")
}

step_signature <- build_step_signature(
  input_files = input_file,
  extra = list(
    step = "05_pseudobulk",
    pseudobulk = cfg$pseudobulk,
    analysis_exclusions = analysis_exclusions
  ),
  code_files = c(
    script_file(),
    file.path(PIPELINE_ROOT, "scripts", "_common.R"),
    file.path(PIPELINE_ROOT, "config", "config.R"),
    analysis_exclusions_file,
    analysis_exclusion_helpers_file
  )
)
initialize_step(
  step_signature,
  checkpoint_dir,
  owned_paths = c(pb_dir, da_dir)
)
if (step_is_complete(
  pb_dir,
  step_signature,
  c(
    all_genes_file,
    significant_file,
    pseudobulk_metadata_file,
    replication_file,
    activity_audit_file,
    cell_exclusion_file,
    cell_exclusion_summary_file,
    abundance_results_file,
    donor_proportions_file,
    abundance_plot_file
  )
)) {
  log_message("Step 05 is already complete and its signature matches. Use --force to rebuild.")
  quit(save = "no", status = 0L)
}
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(da_dir, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(input_file)
counts <- get_assay_counts(object)
metadata <- object@meta.data
metadata$cell_id <- rownames(metadata)

exclusion_result <- mask_excluded_cluster_labels(
  metadata = metadata,
  excluded_clusters = analysis_exclusions$seurat_clusters,
  label_columns = c("cell_type_major", "cell_type_analysis"),
  exclusion_reason = analysis_exclusions$reason
)
metadata <- exclusion_result$metadata
excluded_metadata <- exclusion_result$excluded_metadata

expected_exclusion_counts <- analysis_exclusions$expected_cluster_counts %||% integer(0)
actual_exclusion_counts <- table(as.character(excluded_metadata$seurat_clusters))
if (length(expected_exclusion_counts) > 0L) {
  expected_names <- names(expected_exclusion_counts)
  if (is.null(expected_names) || any(!nzchar(expected_names))) {
    stop("analysis_exclusions$expected_cluster_counts must be a named vector.")
  }
  observed <- as.integer(actual_exclusion_counts[expected_names])
  observed[is.na(observed)] <- 0L
  mismatch <- observed != as.integer(expected_exclusion_counts)
  if (any(mismatch)) {
    stop(
      "Frozen cluster exclusion count mismatch: ",
      paste0(
        expected_names[mismatch], " expected=", as.integer(expected_exclusion_counts[mismatch]),
        " observed=", observed[mismatch], collapse = "; "
      ),
      ". Review the annotated object before running step 05."
    )
  }
}

exclusion_columns <- intersect(
  c(
    "cell_id", "seurat_clusters", "sample_id", "subject_id", "cohort", "condition",
    "cell_type_major", "cell_type_analysis", "analysis_exclusion_reason"
  ),
  colnames(excluded_metadata)
)
safe_write_tsv(excluded_metadata[, exclusion_columns, drop = FALSE], cell_exclusion_file)
exclusion_summary <- data.frame(
  seurat_clusters = names(actual_exclusion_counts),
  n_cells = as.integer(actual_exclusion_counts),
  expected_n_cells = as.integer(expected_exclusion_counts[names(actual_exclusion_counts)]),
  exclusion_reason = analysis_exclusions$reason,
  stringsAsFactors = FALSE
)
safe_write_tsv(exclusion_summary, cell_exclusion_summary_file)
log_message(
  "Excluded ", sum(exclusion_result$excluded), " cells from primary step 05 analyses: cluster(s) ",
  paste(analysis_exclusions$seurat_clusters, collapse = ", ")
)

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) NA_character_ else x[[1]]
}

sample_metadata <- do.call(rbind, lapply(split(metadata, metadata$sample_id), function(x) {
  data.frame(
    sample_id = first_non_missing(x$sample_id),
    subject_id = first_non_missing(x$subject_id),
    cohort = first_non_missing(x$cohort),
    condition = first_non_missing(x$condition),
    age = numeric_or_na(first_non_missing(x$age)),
    sex = first_non_missing(x$sex),
    batch = first_non_missing(x$batch),
    SLEDAI = numeric_or_na(first_non_missing(x$SLEDAI)),
    published_activity_group = first_non_missing(x$sle_activity),
    stringsAsFactors = FALSE
  )
}))
rownames(sample_metadata) <- NULL
sample_metadata$analysis_activity_group <- derive_sledai_group(sample_metadata$SLEDAI)
sample_metadata$activity_binary <- ifelse(
  sample_metadata$analysis_activity_group %in% c("Low", "High"),
  sample_metadata$analysis_activity_group,
  NA_character_
)

published_normalized <- ifelse(
  sample_metadata$published_activity_group == "0 to 3", "Low",
  ifelse(sample_metadata$published_activity_group == ">6", "High",
         ifelse(sample_metadata$published_activity_group == "4 to 6", "Intermediate", NA_character_))
)
activity_audit <- sample_metadata
activity_audit$published_group_normalized <- published_normalized
activity_audit$group_concordant <- is.na(published_normalized) |
  is.na(activity_audit$analysis_activity_group) |
  published_normalized == activity_audit$analysis_activity_group
safe_write_tsv(activity_audit, activity_audit_file)
if (any(!activity_audit$group_concordant)) {
  warning(
    "Published activity strata disagree with numeric SLEDAI for: ",
    paste(activity_audit$sample_id[!activity_audit$group_concordant], collapse = ", "),
    ". Numeric SLEDAI-derived groups are used for analysis."
  )
}

coerce_covariate <- function(x, name) {
  if (name %in% c("age", "SLEDAI")) return(numeric_or_na(x))
  factor(as.character(x))
}

build_adjusted_design <- function(
    data,
    effect_column,
    effect_type = c("binary", "continuous"),
    reference_group = NULL,
    test_group = NULL,
    covariates = character(0)) {
  effect_type <- match.arg(effect_type)
  model_data <- data
  model_data$.row_index <- seq_len(nrow(model_data))
  if (effect_type == "binary") {
    model_data$.effect <- factor(
      as.character(model_data[[effect_column]]),
      levels = c(reference_group, test_group)
    )
  } else {
    model_data$.effect <- numeric_or_na(model_data[[effect_column]])
  }
  model_data <- model_data[!is.na(model_data$.effect), , drop = FALSE]
  if (nrow(model_data) == 0L) return(NULL)

  used <- character(0)
  current_data <- model_data
  current_design <- stats::model.matrix(~.effect, data = current_data)
  if (qr(current_design)$rank < ncol(current_design)) return(NULL)

  for (covariate in covariates) {
    if (!covariate %in% colnames(model_data)) next
    candidate <- model_data
    candidate[[covariate]] <- coerce_covariate(candidate[[covariate]], covariate)
    candidate_used <- c(used, covariate)
    complete <- stats::complete.cases(candidate[c(candidate_used, ".effect")])
    candidate <- candidate[complete, , drop = FALSE]
    if (nrow(candidate) == 0L) next
    for (name in candidate_used) {
      candidate[[name]] <- coerce_covariate(candidate[[name]], name)
    }
    if (any(vapply(candidate[candidate_used], function(x) length(unique(x)) < 2L, logical(1)))) {
      next
    }
    formula <- stats::reformulate(c(candidate_used, ".effect"))
    design <- tryCatch(stats::model.matrix(formula, data = candidate), error = function(e) NULL)
    if (is.null(design) || qr(design)$rank < ncol(design)) next
    used <- candidate_used
    current_data <- candidate
    current_design <- design
  }

  effect_columns <- grep("^\\.effect", colnames(current_design))
  if (length(effect_columns) != 1L) return(NULL)
  list(
    data = current_data,
    design = current_design,
    coefficient = effect_columns[[1]],
    covariates_used = used,
    formula = paste("~", paste(c(used, effect_column), collapse = " + "))
  )
}

aggregate_pseudobulk <- function(
    count_matrix,
    cell_metadata,
    label_column,
    min_cells) {
  labels <- as.character(cell_metadata[[label_column]])
  valid <- !is.na(labels) & nzchar(labels) & !grepl("Unresolved", labels, ignore.case = TRUE)
  valid[is.na(valid)] <- FALSE
  if (!any(valid)) stop("No cells available for pseudobulk level: ", label_column)

  md <- cell_metadata[valid, , drop = FALSE]
  md$cell_label <- as.character(md[[label_column]])
  cell_keys <- paste(md$sample_id, md$cell_label, sep = "||")
  group_frame <- unique(data.frame(
    group_key = cell_keys,
    sample_id = md$sample_id,
    cell_label = md$cell_label,
    stringsAsFactors = FALSE
  ))
  group_frame <- group_frame[order(group_frame$sample_id, group_frame$cell_label), ]
  group_index <- match(cell_keys, group_frame$group_key)
  cell_to_group <- Matrix::sparseMatrix(
    i = seq_along(group_index), j = group_index, x = 1,
    dims = c(length(group_index), nrow(group_frame))
  )
  pb_counts <- count_matrix[, valid, drop = FALSE] %*% cell_to_group
  pb_counts <- methods::as(pb_counts, "dgCMatrix")
  group_frame$n_cells <- tabulate(group_index, nbins = nrow(group_frame))
  keep <- group_frame$n_cells >= min_cells
  group_frame <- group_frame[keep, , drop = FALSE]
  pb_counts <- pb_counts[, keep, drop = FALSE]
  colnames(pb_counts) <- make.unique(
    paste(group_frame$sample_id, make_safe_name(group_frame$cell_label), sep = "__")
  )
  group_frame$pseudobulk_id <- colnames(pb_counts)
  group_frame <- merge(group_frame, sample_metadata, by = "sample_id", all.x = TRUE)
  group_frame <- group_frame[match(colnames(pb_counts), group_frame$pseudobulk_id), , drop = FALSE]
  list(counts = pb_counts, metadata = group_frame)
}

run_edger_binary <- function(
    pb_counts,
    pb_metadata,
    comparison_id,
    group_column,
    reference_group,
    test_group,
    covariates,
    cohort = NULL,
    require_condition = NULL) {
  md <- pb_metadata
  keep <- md[[group_column]] %in% c(reference_group, test_group)
  if (!is.null(cohort)) keep <- keep & md$cohort == cohort
  if (!is.null(require_condition)) keep <- keep & md$condition == require_condition
  keep[is.na(keep)] <- FALSE
  md <- md[keep, , drop = FALSE]
  matrix_subset <- pb_counts[, keep, drop = FALSE]
  if (nrow(md) == 0L) return(data.frame())

  result_list <- list()
  for (label in sort(unique(md$cell_label))) {
    label_keep <- md$cell_label == label
    label_md <- md[label_keep, , drop = FALSE]
    label_counts <- matrix_subset[, label_keep, drop = FALSE]
    design_info <- build_adjusted_design(
      label_md, group_column, "binary", reference_group, test_group, covariates
    )
    if (is.null(design_info)) next
    selected <- design_info$data$.row_index
    label_md <- label_md[selected, , drop = FALSE]
    label_counts <- label_counts[, selected, drop = FALSE]
    group <- factor(label_md[[group_column]], levels = c(reference_group, test_group))
    group_counts <- table(group)
    if (any(group_counts < cfg$pseudobulk$min_samples_per_group)) {
      warning("Skip ", comparison_id, " / ", label, ": insufficient donors after design filtering.")
      next
    }

    table_result <- tryCatch({
      y <- edgeR::DGEList(counts = label_counts, group = group)
      keep_genes <- edgeR::filterByExpr(y, design = design_info$design)
      y <- y[keep_genes, , keep.lib.sizes = FALSE]
      if (nrow(y) < cfg$pseudobulk$min_expressed_genes) {
        stop("too few expressed genes after filterByExpr")
      }
      y <- edgeR::calcNormFactors(y, method = "TMM")
      y <- edgeR::estimateDisp(y, design = design_info$design, robust = TRUE)
      fit <- edgeR::glmQLFit(y, design = design_info$design, robust = TRUE)
      test <- edgeR::glmQLFTest(fit, coef = design_info$coefficient)
      edgeR::topTags(test, n = Inf, sort.by = "PValue")$table
    }, error = function(e) {
      warning("Skip ", comparison_id, " / ", label, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(table_result)) next
    table_result$gene <- rownames(table_result)
    rownames(table_result) <- NULL
    table_result$comparison <- comparison_id
    table_result$cell_label <- label
    table_result$reference_group <- reference_group
    table_result$test_group <- test_group
    table_result$n_reference <- as.integer(group_counts[[reference_group]])
    table_result$n_test <- as.integer(group_counts[[test_group]])
    table_result$covariates_used <- paste(design_info$covariates_used, collapse = ";")
    table_result$design_formula <- design_info$formula
    table_result$significant <- table_result$FDR < cfg$pseudobulk$fdr &
      abs(table_result$logFC) >= cfg$pseudobulk$logfc
    result_list[[label]] <- table_result
  }
  if (length(result_list) == 0L) return(data.frame())
  combine_result_frames(result_list)
}

run_edger_continuous <- function(
    pb_counts,
    pb_metadata,
    comparison_id,
    effect_column,
    covariates,
    cohort,
    require_condition) {
  keep <- pb_metadata$cohort == cohort & pb_metadata$condition == require_condition
  keep <- keep & !is.na(pb_metadata[[effect_column]])
  keep[is.na(keep)] <- FALSE
  md <- pb_metadata[keep, , drop = FALSE]
  matrix_subset <- pb_counts[, keep, drop = FALSE]
  results <- list()
  for (label in sort(unique(md$cell_label))) {
    label_keep <- md$cell_label == label
    label_md <- md[label_keep, , drop = FALSE]
    label_counts <- matrix_subset[, label_keep, drop = FALSE]
    design_info <- build_adjusted_design(
      label_md, effect_column, "continuous", covariates = covariates
    )
    if (is.null(design_info) || nrow(design_info$data) < cfg$pseudobulk$min_samples_continuous) next
    selected <- design_info$data$.row_index
    label_counts <- label_counts[, selected, drop = FALSE]
    table_result <- tryCatch({
      y <- edgeR::DGEList(counts = label_counts)
      keep_genes <- edgeR::filterByExpr(y, design = design_info$design)
      y <- y[keep_genes, , keep.lib.sizes = FALSE]
      if (nrow(y) < cfg$pseudobulk$min_expressed_genes) {
        stop("too few expressed genes after filterByExpr")
      }
      y <- edgeR::calcNormFactors(y, method = "TMM")
      y <- edgeR::estimateDisp(y, design = design_info$design, robust = TRUE)
      fit <- edgeR::glmQLFit(y, design = design_info$design, robust = TRUE)
      test <- edgeR::glmQLFTest(fit, coef = design_info$coefficient)
      edgeR::topTags(test, n = Inf, sort.by = "PValue")$table
    }, error = function(e) {
      warning("Skip ", comparison_id, " / ", label, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(table_result)) next
    table_result$gene <- rownames(table_result)
    rownames(table_result) <- NULL
    table_result$comparison <- comparison_id
    table_result$cell_label <- label
    table_result$effect <- effect_column
    table_result$effect_scale <- "per one SLEDAI point"
    table_result$n_samples <- nrow(design_info$data)
    table_result$covariates_used <- paste(design_info$covariates_used, collapse = ";")
    table_result$design_formula <- design_info$formula
    table_result$significant <- table_result$FDR < cfg$pseudobulk$fdr &
      abs(table_result$logFC) >= cfg$pseudobulk$logfc
    results[[label]] <- table_result
  }
  if (length(results) == 0L) return(data.frame())
  combine_result_frames(results)
}

build_abundance_table <- function(cell_metadata, label_column) {
  labels <- as.character(cell_metadata[[label_column]])
  valid <- !is.na(labels) & nzchar(labels) & !grepl("Unresolved", labels, ignore.case = TRUE)
  valid[is.na(valid)] <- FALSE
  md <- cell_metadata[valid, , drop = FALSE]
  md$cell_label <- as.character(md[[label_column]])
  observed <- as.data.frame(table(sample_id = md$sample_id, cell_label = md$cell_label), stringsAsFactors = FALSE)
  totals <- as.data.frame(table(sample_id = md$sample_id), stringsAsFactors = FALSE)
  colnames(totals)[[2]] <- "total_cells"
  observed <- merge(observed, totals, by = "sample_id", all.x = TRUE)
  observed <- merge(observed, sample_metadata, by = "sample_id", all.x = TRUE)
  observed$proportion <- observed$Freq / observed$total_cells
  observed
}

fit_equal_weight_linear <- function(outcome, design, coefficient) {
  xtx_inverse <- tryCatch(solve(crossprod(design)), error = function(e) NULL)
  if (is.null(xtx_inverse)) return(NULL)
  beta <- as.numeric(xtx_inverse %*% crossprod(design, outcome))
  residual <- outcome - as.numeric(design %*% beta)
  df <- nrow(design) - ncol(design)
  if (df <= 0L) return(NULL)
  sigma2 <- sum(residual^2) / df
  se <- sqrt(diag(xtx_inverse) * sigma2)
  statistic <- beta[[coefficient]] / se[[coefficient]]
  p_value <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  critical <- stats::qt((1 + cfg$inference$confidence_level) / 2, df = df)
  list(
    estimate = beta[[coefficient]],
    standard_error = se[[coefficient]],
    statistic = statistic,
    degrees_freedom = df,
    p_value = p_value,
    conf_low = beta[[coefficient]] - critical * se[[coefficient]],
    conf_high = beta[[coefficient]] + critical * se[[coefficient]]
  )
}

run_abundance_model <- function(
    abundance,
    comparison_id,
    effect_column,
    effect_type,
    covariates,
    reference_group = NULL,
    test_group = NULL,
    cohort,
    require_condition = NULL) {
  data <- abundance[abundance$cohort == cohort, , drop = FALSE]
  if (!is.null(require_condition)) data <- data[data$condition == require_condition, , drop = FALSE]
  if (effect_type == "binary") {
    data <- data[data[[effect_column]] %in% c(reference_group, test_group), , drop = FALSE]
  } else {
    data <- data[!is.na(data[[effect_column]]), , drop = FALSE]
  }
  results <- lapply(sort(unique(data$cell_label)), function(label) {
    x <- data[data$cell_label == label, , drop = FALSE]
    design_info <- build_adjusted_design(
      x, effect_column, effect_type, reference_group, test_group, covariates
    )
    if (is.null(design_info)) return(NULL)
    if (effect_type == "binary") {
      group_counts <- table(factor(
        design_info$data[[effect_column]], levels = c(reference_group, test_group)
      ))
      if (any(group_counts < cfg$pseudobulk$min_samples_per_group)) return(NULL)
    } else if (nrow(design_info$data) < cfg$pseudobulk$min_samples_continuous) {
      return(NULL)
    }
    selected <- design_info$data$.row_index
    x <- x[selected, , drop = FALSE]
    pseudocount <- cfg$pseudobulk$abundance_logit_pseudocount
    outcome <- log((x$Freq + pseudocount) / (x$total_cells - x$Freq + pseudocount))
    fit <- fit_equal_weight_linear(outcome, design_info$design, design_info$coefficient)
    if (is.null(fit)) return(NULL)
    data.frame(
      comparison = comparison_id,
      cell_label = label,
      effect = effect_column,
      estimate_logit_scale = fit$estimate,
      standard_error = fit$standard_error,
      conf_low = fit$conf_low,
      conf_high = fit$conf_high,
      statistic = fit$statistic,
      degrees_freedom = fit$degrees_freedom,
      p_value = fit$p_value,
      n_samples = nrow(x),
      mean_reference = if (effect_type == "binary") mean(x$proportion[x[[effect_column]] == reference_group]) else NA_real_,
      mean_test = if (effect_type == "binary") mean(x$proportion[x[[effect_column]] == test_group]) else NA_real_,
      covariates_used = paste(design_info$covariates_used, collapse = ";"),
      design_formula = design_info$formula,
      donor_weighting = "equal",
      stringsAsFactors = FALSE
    )
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) return(data.frame())
  results <- do.call(rbind, results)
  results$fdr <- stats::p.adjust(results$p_value, method = "BH")
  results
}

analysis_levels <- list(
  major = list(label_column = "cell_type_major", min_cells = cfg$pseudobulk$min_cells_major),
  analysis = list(label_column = "cell_type_analysis", min_cells = cfg$pseudobulk$min_cells_analysis)
)
if (!all(vapply(analysis_levels, function(x) x$label_column %in% colnames(metadata), logical(1)))) {
  stop("Annotated object is missing one or more configured cell-label columns.")
}

aggregation_cache <- load_checkpoint_rds(aggregation_cache_path, step_signature)
if (is.null(aggregation_cache) || !is.list(aggregation_cache)) aggregation_cache <- list()
for (level_name in names(analysis_levels)) {
  cached_level <- aggregation_cache[[level_name]]
  if (is.list(cached_level) && all(c("counts", "metadata") %in% names(cached_level))) {
    log_message("Reusing pseudobulk aggregation level: ", level_name)
    next
  }
  level <- analysis_levels[[level_name]]
  log_message("Pseudobulk aggregation level: ", level_name)
  aggregation_cache[[level_name]] <- aggregate_pseudobulk(
    counts,
    metadata,
    level$label_column,
    level$min_cells
  )
  save_checkpoint_rds(
    aggregation_cache,
    aggregation_cache_path,
    step_signature,
    compress = TRUE
  )
}
if (!named_cache_has_fields(
  aggregation_cache,
  names(analysis_levels),
  c("counts", "metadata")
)) {
  stop("Pseudobulk aggregation cache is incomplete after aggregation.")
}

all_de_results <- list()
all_da_results <- list()
all_pb_metadata <- list()
all_replication <- list()
all_abundance <- list()
for (level_name in names(analysis_levels)) {
  level <- analysis_levels[[level_name]]
  pb <- aggregation_cache[[level_name]]
  level_metadata <- pb$metadata
  level_metadata$analysis_level <- level_name
  all_pb_metadata[[level_name]] <- level_metadata

  child_result <- run_edger_binary(
    pb$counts, pb$metadata, "Child_SLE_vs_HD", "condition", "HD", "SLE",
    cfg$pseudobulk$covariates_child, cohort = "Child"
  )
  adult_result <- run_edger_binary(
    pb$counts, pb$metadata, "Adult_SLE_vs_HD", "condition", "HD", "SLE",
    cfg$pseudobulk$covariates_adult, cohort = "Adult"
  )
  activity_result <- run_edger_binary(
    pb$counts, pb$metadata, "Child_SLE_High_vs_Low_Activity", "activity_binary", "Low", "High",
    cfg$pseudobulk$covariates_activity, cohort = "Child", require_condition = "SLE"
  )
  sledai_result <- run_edger_continuous(
    pb$counts, pb$metadata, "Child_SLE_SLEDAI_continuous", "SLEDAI",
    cfg$pseudobulk$covariates_activity, "Child", "SLE"
  )
  level_results <- list(child_result, adult_result, activity_result, sledai_result)
  level_results <- lapply(level_results, function(result) {
    if (is.data.frame(result) && nrow(result) > 0L) result$analysis_level <- level_name
    result
  })
  child_result <- level_results[[1]]
  adult_result <- level_results[[2]]
  activity_result <- level_results[[3]]
  sledai_result <- level_results[[4]]
  all_de_results[[paste0(level_name, "_child")]] <- child_result
  all_de_results[[paste0(level_name, "_adult")]] <- adult_result
  all_de_results[[paste0(level_name, "_activity")]] <- activity_result
  all_de_results[[paste0(level_name, "_sledai")]] <- sledai_result

  if (nrow(child_result) > 0L && nrow(adult_result) > 0L) {
    replication <- merge(
      child_result[c("gene", "cell_label", "logFC", "FDR")],
      adult_result[c("gene", "cell_label", "logFC", "FDR")],
      by = c("gene", "cell_label"), suffixes = c("_child", "_adult")
    )
    replication$direction_concordant <- sign(replication$logFC_child) == sign(replication$logFC_adult)
    replication$child_significant <- replication$FDR_child < cfg$pseudobulk$fdr &
      abs(replication$logFC_child) >= cfg$pseudobulk$logfc
    replication$adult_support <- replication$FDR_adult < cfg$pseudobulk$adult_support_fdr
    replication$analysis_level <- level_name
    all_replication[[level_name]] <- replication
  }

  abundance <- build_abundance_table(metadata, level$label_column)
  abundance$analysis_level <- level_name
  all_abundance[[level_name]] <- abundance
  da_results <- list(
    run_abundance_model(
      abundance, "Child_SLE_vs_HD", "condition", "binary",
      cfg$pseudobulk$covariates_child, "HD", "SLE", "Child"
    ),
    run_abundance_model(
      abundance, "Adult_SLE_vs_HD", "condition", "binary",
      cfg$pseudobulk$covariates_adult, "HD", "SLE", "Adult"
    ),
    run_abundance_model(
      abundance, "Child_SLE_High_vs_Low_Activity", "activity_binary", "binary",
      cfg$pseudobulk$covariates_activity, "Low", "High", "Child", "SLE"
    ),
    run_abundance_model(
      abundance, "Child_SLE_SLEDAI_continuous", "SLEDAI", "continuous",
      cfg$pseudobulk$covariates_activity, cohort = "Child", require_condition = "SLE"
    )
  )
  da_combined <- combine_result_frames(da_results)
  if (nrow(da_combined) > 0L) {
    da_combined$analysis_level <- level_name
    all_da_results[[level_name]] <- da_combined
  }

  if (level_name == "analysis") {
    plot_data <- abundance[abundance$cohort == "Child", , drop = FALSE]
    proportion_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = cell_label, y = proportion, fill = condition)
    ) +
      ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge(width = 0.75)) +
      ggplot2::geom_point(
        ggplot2::aes(color = condition),
        position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
        alpha = 0.7, size = 1
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "Child cohort: donor-level cell proportions", x = NULL, y = "Cell fraction")
    ggplot2::ggsave(
      abundance_plot_file,
      proportion_plot, width = 12, height = 6
    )
  }
}

all_de_combined <- combine_result_frames(all_de_results)
if (nrow(all_de_combined) == 0L) {
  all_de_combined <- data.frame(
    status = "no_estimable_pseudobulk_results",
    analysis_level = NA_character_,
    stringsAsFactors = FALSE
  )
}
atomic_save_rds(all_de_combined, all_genes_file, compress = TRUE)
if ("significant" %in% colnames(all_de_combined)) {
  significant <- all_de_combined[
    !is.na(all_de_combined$significant) & all_de_combined$significant,
    ,
    drop = FALSE
  ]
} else {
  significant <- data.frame(
    status = "no_estimable_pseudobulk_results",
    analysis_level = NA_character_,
    stringsAsFactors = FALSE
  )
}
safe_write_tsv(significant, significant_file)
safe_write_tsv(combine_result_frames(all_pb_metadata), pseudobulk_metadata_file)

replication_combined <- combine_result_frames(all_replication)
if (nrow(replication_combined) == 0L) {
  replication_combined <- data.frame(
    status = "no_child_adult_replication_results",
    analysis_level = NA_character_,
    stringsAsFactors = FALSE
  )
}
safe_write_tsv(replication_combined, replication_file)

abundance_combined <- combine_result_frames(all_abundance)
safe_write_tsv(abundance_combined, donor_proportions_file)
da_combined <- combine_result_frames(all_da_results)
if (nrow(da_combined) == 0L) {
  da_combined <- data.frame(
    status = "no_estimable_abundance_results",
    analysis_level = NA_character_,
    stringsAsFactors = FALSE
  )
}
safe_write_tsv(da_combined, abundance_results_file)

if (!file.exists(abundance_plot_file)) {
  placeholder <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No estimable abundance plot")
  ggplot2::ggsave(abundance_plot_file, placeholder, width = 8, height = 4)
}

write_session_info(file.path(cfg$paths$logs, "05_pseudobulk_sessionInfo.txt"))
critical_outputs <- c(
  aggregation_cache_path,
  all_genes_file,
  significant_file,
  pseudobulk_metadata_file,
  replication_file,
  activity_audit_file,
  cell_exclusion_file,
  cell_exclusion_summary_file,
  abundance_results_file,
  donor_proportions_file,
  abundance_plot_file
)
legacy_outputs <- c(
  file.path(pb_dir, c("major", "analysis", "author_subcluster")),
  file.path(da_dir, c("major", "analysis", "author_subcluster")),
  file.path(pb_dir, "all_pseudobulk_results.tsv"),
  file.path(da_dir, "all_abundance_results.tsv")
)
complete_step_and_cleanup(
  pb_dir,
  step_signature,
  critical_outputs,
  cleanup_paths = legacy_outputs
)
log_message("Pseudobulk and differential abundance analyses completed.")
