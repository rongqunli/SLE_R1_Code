if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(c("data.table", "limma", "matrixStats", "WGCNA"),
                 "longitudinal WGCNA")
suppressPackageStartupMessages(library(WGCNA))

script_path <- script_file()
expression_file <- file.path(cfg$paths$processed, "01_gene_expression.rds")
clinical_file <- file.path(
  cfg$paths$results, "03_longitudinal", "longitudinal_analysis_data.tsv"
)
output_dir <- file.path(cfg$paths$results, "06_wgcna")
object_file <- file.path(cfg$paths$processed, "04_wgcna.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- c(
  object_file,
  file.path(output_dir, c(
    "module_assignments.tsv", "module_eigengenes.tsv",
    "module_trait_associations.tsv", "module_membership.tsv",
    "hub_genes.tsv", "soft_threshold_diagnostics.tsv",
    "module_stability.tsv", "wgcna_input_audit.tsv"
  ))
)
signature <- build_step_signature(
  input_files = c(expression_file, clinical_file),
  code_files = c(script_path, file.path(PIPELINE_ROOT, "scripts", "_common.R"),
                 file.path(PIPELINE_ROOT, "config", "config.R")),
  extra = list(
    method = "patient_centered_signed_WGCNA_v1", selected_genes = 5000L,
    bootstrap_iterations = 25L
  )
)
if (!has_flag("--force") && step_is_complete(output_dir, signature, outputs)) {
  log_message("Step 06 WGCNA is already complete.")
  quit(save = "no", status = 0L)
}

center_expression_by_patient <- function(expression, patient_id) {
  centered <- expression
  patient_columns <- split(seq_along(patient_id), patient_id)
  for (columns in patient_columns) {
    centered[, columns] <- sweep(
      expression[, columns, drop = FALSE], 1L,
      rowMeans(expression[, columns, drop = FALSE]), FUN = "-"
    )
  }
  centered
}

fit_module_family <- function(
    eigengenes, data, axis, family, rhs, target_term, predictor) {
  design_formula <- stats::as.formula(paste("~", rhs))
  variables <- all.vars(design_formula)
  complete <- stats::complete.cases(data[, ..variables]) & !is.na(data$patient_id)
  model_data <- droplevels(as.data.frame(data[complete]))
  model_eigengenes <- eigengenes[complete, , drop = FALSE]
  design <- stats::model.matrix(design_formula, data = model_data)
  if (qr(design)$rank != ncol(design)) {
    stop("Rank-deficient module design for model family: ", family)
  }
  if (!target_term %in% colnames(design)) {
    stop("Target term missing from module design: ", target_term)
  }
  block <- factor(model_data$patient_id)
  correlation <- tryCatch(
    limma::duplicateCorrelation(
      t(model_eigengenes), design = design, block = block
    )$consensus.correlation,
    error = function(error) {
      log_message("duplicateCorrelation failed for ", family, ": ",
                  conditionMessage(error), "; using zero correlation.")
      0
    }
  )
  if (!is.finite(correlation)) correlation <- 0
  fit <- limma::lmFit(
    t(model_eigengenes), design = design, block = block,
    correlation = correlation
  )
  fit <- limma::eBayes(fit, trend = TRUE)
  target <- match(target_term, colnames(design))
  predictor_sd <- stats::sd(as.numeric(predictor[complete]), na.rm = TRUE)
  if (!is.finite(predictor_sd) || predictor_sd <= 0) {
    stop("Module predictor has no finite variation for: ", family)
  }
  standard_error <- sqrt(fit$s2.post) * fit$stdev.unscaled[, target]
  degrees_freedom <- fit$df.total
  critical_value <- stats::qt(0.975, df = degrees_freedom)
  estimates <- unname(fit$coefficients[, target])
  p_values <- unname(fit$p.value[, target])
  data.table::data.table(
    axis = axis,
    model_family = family,
    module = sub("^ME", "", rownames(fit$coefficients)),
    term = target_term,
    estimate = estimates,
    std_error = standard_error,
    conf_low = estimates - critical_value * standard_error,
    conf_high = estimates + critical_value * standard_error,
    t_statistic = unname(fit$t[, target]),
    p_value = p_values,
    fdr = stats::p.adjust(p_values, method = "BH"),
    effect_per_sd = estimates * predictor_sd,
    predictor_sd = predictor_sd,
    n_patients = data.table::uniqueN(model_data$patient_id),
    n_visits = nrow(model_data),
    consensus_correlation = correlation,
    formula = paste0(
      "module eigengene ~ ", rhs,
      " + patient block (limma duplicateCorrelation)"
    )
  )
}

bootstrap_membership_stability <- function(
    expression, eigengenes, modules, patients, n_bootstrap, seed) {
  module_names <- colnames(eigengenes)
  module_colors <- sub("^ME", "", module_names)
  patient_rows <- split(seq_along(patients), patients)
  patient_levels <- names(patient_rows)
  full_membership <- stats::cor(
    expression, eigengenes, use = "pairwise.complete.obs"
  )
  correlations <- matrix(
    NA_real_, nrow = n_bootstrap, ncol = length(module_colors),
    dimnames = list(NULL, module_colors)
  )
  set.seed(seed)
  for (iteration in seq_len(n_bootstrap)) {
    sampled_patients <- sample(
      patient_levels, length(patient_levels), replace = TRUE
    )
    rows <- unlist(patient_rows[sampled_patients], use.names = FALSE)
    for (module_index in seq_along(module_colors)) {
      genes <- which(modules == module_colors[[module_index]])
      bootstrap_membership <- stats::cor(
        expression[rows, genes, drop = FALSE],
        eigengenes[rows, module_index],
        use = "pairwise.complete.obs"
      )[, 1L]
      correlations[iteration, module_index] <- suppressWarnings(stats::cor(
        full_membership[genes, module_index], bootstrap_membership,
        method = "spearman", use = "complete.obs"
      ))
    }
  }
  data.table::rbindlist(lapply(seq_along(module_colors), function(index) {
    values <- correlations[, index]
    values <- values[is.finite(values)]
    data.table::data.table(
      module = module_colors[[index]],
      n_genes = sum(modules == module_colors[[index]]),
      n_bootstrap = n_bootstrap,
      n_finite_bootstrap = length(values),
      median_kme_rank_correlation = stats::median(values),
      q025_kme_rank_correlation = unname(stats::quantile(values, 0.025)),
      q975_kme_rank_correlation = unname(stats::quantile(values, 0.975)),
      proportion_correlation_ge_0_7 = mean(values >= 0.7),
      stability_method = paste(
        "patient bootstrap; fixed modules; rank correlation of gene kME"
      )
    )
  }), use.names = TRUE)
}

expression_object <- readRDS(expression_file)
expression <- expression_object$expression
clinical <- data.table::fread(clinical_file, sep = "\t", na.strings = "NA")
if (!is.matrix(expression) || nrow(expression) != 19307L) {
  stop("Unexpected expression matrix for WGCNA.")
}
repeated_data <- data.table::copy(clinical[patient_visit_count >= 2L])
if (nrow(repeated_data) != 897L || data.table::uniqueN(repeated_data$patient_id) != 131L) {
  stop("Unexpected repeated-visit WGCNA cohort.")
}
columns <- match(repeated_data$sample_id, colnames(expression))
if (anyNA(columns)) stop("WGCNA samples are missing from expression matrix.")
repeated_expression <- expression[, columns, drop = FALSE]
repeated_data[, `:=`(
  patient_id = factor(patient_id),
  batch = factor(as.character(batch)),
  renal_numeric = as.numeric(renal_activity_numeric)
)]
repeated_data[, renal_between := mean(renal_numeric, na.rm = TRUE), by = patient_id]
repeated_data[, renal_within := renal_numeric - renal_between]

log_message("Patient-centering and batch-residualizing 19,307 genes.")
centered <- center_expression_by_patient(
  repeated_expression, repeated_data$patient_id
)
batch_residualized <- limma::removeBatchEffect(
  centered, batch = repeated_data$batch
)
batch_residualized <- center_expression_by_patient(
  batch_residualized, repeated_data$patient_id
)
selection_variance <- matrixStats::rowVars(batch_residualized)
selected_indices <- head(order(selection_variance, decreasing = TRUE), 5000L)
selected_expression <- batch_residualized[selected_indices, , drop = FALSE]
dat_expr <- t(selected_expression)
rownames(dat_expr) <- repeated_data$sample_id

quality <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
if (!quality$allOK) {
  stop("Selected WGCNA matrix contains unsuitable genes or samples.")
}

powers <- 1:12
log_message("Evaluating signed-network soft-threshold powers 1-12.")
soft_threshold <- WGCNA::pickSoftThreshold(
  dat_expr, powerVector = powers, networkType = "signed",
  corFnc = WGCNA::bicor, corOptions = list(maxPOutliers = 0.1), verbose = 3
)
fit_indices <- data.table::as.data.table(soft_threshold$fitIndices)
data.table::setnames(fit_indices, make.names(names(fit_indices), unique = TRUE))
power_column <- names(fit_indices)[tolower(names(fit_indices)) == "power"][[1]]
r2_column <- grep("SFT.*R.*sq", names(fit_indices), value = TRUE)[[1]]
slope_column <- grep("slope", names(fit_indices), ignore.case = TRUE, value = TRUE)[[1]]
eligible <- which(fit_indices[[r2_column]] >= 0.8 & fit_indices[[slope_column]] < 0)
if (length(eligible)) {
  chosen_row <- eligible[[1]]
} else {
  finite_rows <- which(is.finite(fit_indices[[r2_column]]))
  chosen_row <- finite_rows[[which.max(fit_indices[[r2_column]][finite_rows])]]
}
chosen_power <- as.integer(fit_indices[[power_column]][chosen_row])
fit_indices[, selected := get(power_column) == chosen_power]
data.table::setnames(fit_indices, power_column, "power")

log_message("Building signed WGCNA network at power ", chosen_power, ".")
set.seed(cfg$seed)
network <- WGCNA::blockwiseModules(
  dat_expr, power = chosen_power, maxBlockSize = 6000L,
  corType = "bicor", maxPOutliers = 0.1,
  networkType = "signed", TOMType = "signed",
  minModuleSize = 30L, deepSplit = 2L, mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE, numericLabels = TRUE,
  saveTOMs = FALSE, randomSeed = cfg$seed, verbose = 3
)
module_colors <- WGCNA::labels2colors(network$colors)
if (length(setdiff(unique(module_colors), "grey")) < 2L) {
  stop("WGCNA found fewer than two non-grey modules.")
}
eigengenes <- WGCNA::moduleEigengenes(
  dat_expr, colors = module_colors, align = "along average"
)$eigengenes
eigengenes <- WGCNA::orderMEs(eigengenes)
eigengenes <- eigengenes[, sub("^ME", "", colnames(eigengenes)) != "grey", drop = FALSE]
module_names <- sub("^ME", "", colnames(eigengenes))

assignments <- data.table::data.table(
  gene = colnames(dat_expr),
  module = module_colors,
  module_label = as.integer(network$colors),
  selection_variance = selection_variance[selected_indices]
)
assignments[, module_size := .N, by = module]

all_kme <- stats::cor(dat_expr, eigengenes, use = "pairwise.complete.obs")
own_column <- match(module_colors, module_names)
own_kme <- rep(0, nrow(assignments))
non_grey_gene <- which(!is.na(own_column))
own_kme[non_grey_gene] <- all_kme[
  cbind(non_grey_gene, own_column[non_grey_gene])
]
best_index <- max.col(abs(all_kme), ties.method = "first")
membership <- data.table::data.table(
  gene = assignments$gene,
  module = assignments$module,
  kME = own_kme,
  abs_kME = abs(own_kme),
  best_correlated_module = module_names[best_index],
  best_module_kME = all_kme[cbind(seq_len(nrow(all_kme)), best_index)],
  selection_variance = assignments$selection_variance
)
membership[, hub_rank := data.table::frank(-abs_kME, ties.method = "first"), by = module]
hubs <- membership[module != "grey" & hub_rank <= 25L]
data.table::setorder(hubs, module, hub_rank, gene)

composition_covariates <- c(
  "NK", "T_CELL", "NEUTROPHIL", "MONOCYTE", "ERYTHROID", "PLATELET", "B_CELL"
)
composition_rhs <- paste(composition_covariates, collapse = " + ")
specifications <- list(
  list(axis = "general_activity", family = "activity_base",
       rhs = "SLEDAI_within + SLEDAI_between + timepoint_scaled + batch",
       term = "SLEDAI_within", predictor = repeated_data$SLEDAI_within),
  list(axis = "general_activity", family = "activity_composition_adjusted",
       rhs = paste("SLEDAI_within + SLEDAI_between + timepoint_scaled + batch +", composition_rhs),
       term = "SLEDAI_within", predictor = repeated_data$SLEDAI_within),
  list(axis = "renal_complement", family = "renal_base",
       rhs = "renal_within + renal_between + timepoint_scaled + batch",
       term = "renal_within", predictor = repeated_data$renal_within),
  list(axis = "renal_complement", family = "renal_composition_adjusted",
       rhs = paste("renal_within + renal_between + timepoint_scaled + batch +", composition_rhs),
       term = "renal_within", predictor = repeated_data$renal_within),
  list(axis = "renal_complement", family = "c3_base",
       rhs = "C3_within + C3_between + timepoint_scaled + batch",
       term = "C3_within", predictor = repeated_data$C3_within),
  list(axis = "renal_complement", family = "c3_composition_adjusted",
       rhs = paste("C3_within + C3_between + timepoint_scaled + batch +", composition_rhs),
       term = "C3_within", predictor = repeated_data$C3_within),
  list(axis = "renal_complement", family = "c4_base",
       rhs = "C4_within + C4_between + timepoint_scaled + batch",
       term = "C4_within", predictor = repeated_data$C4_within),
  list(axis = "renal_complement", family = "c4_composition_adjusted",
       rhs = paste("C4_within + C4_between + timepoint_scaled + batch +", composition_rhs),
       term = "C4_within", predictor = repeated_data$C4_within)
)

association_list <- lapply(specifications, function(specification) {
  log_message("Fitting module model family ", specification$family, ".")
  fit_module_family(
    eigengenes, repeated_data, specification$axis, specification$family,
    specification$rhs, specification$term, specification$predictor
  )
})
associations <- data.table::rbindlist(association_list, use.names = TRUE)

log_message("Estimating module-membership stability with 25 patient bootstraps.")
stability <- bootstrap_membership_stability(
  dat_expr, eigengenes, module_colors, repeated_data$patient_id,
  n_bootstrap = 25L, seed = cfg$seed + 6L
)

eigengene_table <- data.table::as.data.table(eigengenes)
eigengene_table <- cbind(
  repeated_data[, .(sample_id, patient_id)], eigengene_table
)
audit <- data.table::data.table(
  metric = c(
    "input_visits", "input_patients", "input_genes", "selected_genes",
    "non_grey_modules", "grey_genes", "chosen_power", "bootstrap_iterations"
  ),
  value = c(
    nrow(repeated_data), data.table::uniqueN(repeated_data$patient_id),
    nrow(expression), nrow(assignments), length(module_names),
    sum(module_colors == "grey"), chosen_power, 25L
  )
)

safe_write_tsv(assignments, outputs[[2]])
safe_write_tsv(eigengene_table, outputs[[3]])
safe_write_tsv(associations, outputs[[4]])
safe_write_tsv(membership, outputs[[5]])
safe_write_tsv(hubs, outputs[[6]])
safe_write_tsv(fit_indices, outputs[[7]])
safe_write_tsv(stability, outputs[[8]])
safe_write_tsv(audit, outputs[[9]])
atomic_save_rds(list(
  module_assignments = as.data.frame(assignments),
  module_eigengenes = as.data.frame(eigengene_table),
  associations = as.data.frame(associations),
  membership = as.data.frame(membership),
  stability = as.data.frame(stability),
  soft_threshold = list(chosen_power = chosen_power,
                        diagnostics = as.data.frame(fit_indices)),
  analysis_data = as.data.frame(repeated_data)
), object_file, compress = TRUE)
write_session_info(file.path(cfg$paths$logs, "06_wgcna_sessionInfo.txt"))
mark_step_complete(output_dir, signature, outputs)
log_message(
  "WGCNA complete: ", length(module_names), " non-grey modules across ",
  nrow(assignments), " genes."
)
