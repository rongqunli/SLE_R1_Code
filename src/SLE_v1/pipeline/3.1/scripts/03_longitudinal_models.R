if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(c("data.table", "lme4"), "patient-level longitudinal models")

script_path <- script_file()
score_file <- file.path(cfg$paths$processed, "02_program_scores.rds")
qc_file <- file.path(cfg$paths$results, "01_expression_qc", "sample_qc_metrics.tsv")
output_dir <- file.path(cfg$paths$results, "03_longitudinal")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- file.path(output_dir, c(
  "longitudinal_analysis_data.tsv", "longitudinal_model_results.tsv",
  "primary_leave_one_patient_out.tsv", "model_diagnostics.tsv",
  "anti_dsDNA_format_audit.tsv", "longitudinal_analysis.rds"
))
signature <- build_step_signature(
  input_files = c(score_file, qc_file, cfg$files$programs, cfg$paths$frozen_meta),
  code_files = c(script_path, file.path(PIPELINE_ROOT, "scripts", "_common.R"),
                 file.path(PIPELINE_ROOT, "config", "config.R"))
)
if (!has_flag("--force") && step_is_complete(output_dir, signature, outputs)) {
  log_message("Step 03 longitudinal models are already complete.")
  quit(save = "no", status = 0L)
}
if (!file.exists(score_file)) stop("Missing step 02 program score object.")
if (!file.exists(qc_file)) stop("Missing step 01 sample QC metrics.")

parse_numeric_field <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "not_available", "not_applicable", "Data Not Available", "unknown")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

center_within_patient <- function(data, value, prefix) {
  value_column <- paste0(prefix, "_numeric")
  between_column <- paste0(prefix, "_between")
  within_column <- paste0(prefix, "_within")
  data[, (value_column) := parse_numeric_field(get(value))]
  data[, (between_column) := {
    observed <- get(value_column)
    if (all(is.na(observed))) NA_real_ else mean(observed, na.rm = TRUE)
  }, by = patient_id]
  data[, (within_column) := get(value_column) - get(between_column)]
  invisible(data)
}

sanitize_message <- function(x) {
  x <- paste(unique(x[nzchar(x)]), collapse = " | ")
  gsub("[\r\n\t]+", " ", x)
}

fit_target_term <- function(data, outcome, rhs, model_family, target_term) {
  formula_text <- paste(outcome, "~", rhs)
  formula <- stats::as.formula(formula_text)
  variables <- all.vars(formula)
  model_data <- droplevels(as.data.frame(data[stats::complete.cases(data[, ..variables])]))
  warning_messages <- character(0)
  error_message <- ""
  fit <- tryCatch(
    withCallingHandlers(
      lme4::lmer(
        formula,
        data = model_data,
        REML = FALSE,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 200000),
          check.conv.singular = "ignore"
        )
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )

  base <- data.table::data.table(
    model_family = model_family,
    program = outcome,
    term = target_term,
    estimate = NA_real_,
    std_error = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    p_value = NA_real_,
    n_patients = data.table::uniqueN(model_data$patient_id),
    n_visits = nrow(model_data),
    formula = formula_text,
    converged = FALSE,
    singular = NA,
    warning = sanitize_message(c(warning_messages, error_message))
  )
  diagnostic <- data.table::copy(base)[, c(
    "estimate", "std_error", "conf_low", "conf_high", "p_value"
  ) := NULL]

  if (is.null(fit)) return(list(result = base, diagnostic = diagnostic, fit = NULL))
  coefficients <- stats::coef(summary(fit))
  convergence_messages <- fit@optinfo$conv$lme4$messages
  fit_converged <- is.null(convergence_messages)
  fit_singular <- lme4::isSingular(fit, tol = 1e-4)
  combined_warning <- sanitize_message(c(
    warning_messages,
    if (is.null(convergence_messages)) character(0) else convergence_messages
  ))
  base[, `:=`(
    converged = fit_converged,
    singular = fit_singular,
    warning = combined_warning
  )]
  diagnostic[, `:=`(
    converged = fit_converged,
    singular = fit_singular,
    warning = combined_warning
  )]

  if (!target_term %in% rownames(coefficients)) {
    base[, warning := sanitize_message(c(warning, "Target term absent from fitted coefficients."))]
    diagnostic[, warning := base$warning]
    return(list(result = base, diagnostic = diagnostic, fit = fit))
  }

  beta_value <- unname(coefficients[target_term, "Estimate"])
  se_value <- unname(coefficients[target_term, "Std. Error"])
  z_value <- beta_value / se_value
  base[, `:=`(
    estimate = beta_value,
    std_error = se_value,
    conf_low = beta_value - 1.96 * se_value,
    conf_high = beta_value + 1.96 * se_value,
    p_value = 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
  )]
  list(result = base, diagnostic = diagnostic, fit = fit)
}

score_object <- readRDS(score_file)
metadata <- data.table::as.data.table(score_object$metadata)
score_wide <- data.table::as.data.table(t(score_object$score_matrix), keep.rownames = "sample_id")
analysis <- merge(metadata, score_wide, by = "sample_id", all.x = TRUE, sort = FALSE)
qc_metrics <- data.table::fread(
  qc_file,
  select = c("sample_id", "disease_state", "outlier_flag")
)
if (anyDuplicated(qc_metrics$sample_id)) stop("Sample QC metrics contain duplicate sample IDs.")
qc_metrics[, qc_outlier_flag := as.logical(outlier_flag)]
if (anyNA(qc_metrics$qc_outlier_flag)) stop("Sample QC outlier flags must be TRUE or FALSE.")
flagged_qc <- qc_metrics[qc_outlier_flag == TRUE]
if (nrow(flagged_qc) != 6L || any(flagged_qc$disease_state != "SLE")) {
  stop("Expected exactly 6 QC-flagged samples, all from SLE visits.")
}
analysis <- merge(
  analysis,
  qc_metrics[, .(sample_id, qc_outlier_flag)],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
analysis <- analysis[match(metadata$sample_id, sample_id)]
analysis <- analysis[disease_state == "SLE"]
if (nrow(analysis) != cfg$expected$sle_visits) stop("Unexpected SLE analysis row count.")
if (anyNA(analysis$qc_outlier_flag) || sum(analysis$qc_outlier_flag) != 6L) {
  stop("Expected exactly 6 QC-flagged SLE samples in the analysis data.")
}

center_within_patient(analysis, "SLEDAI", "SLEDAI")
center_within_patient(analysis, "C3", "C3")
center_within_patient(analysis, "C4", "C4")
analysis[, `:=`(
  timepoint_numeric = parse_numeric_field(timepoint),
  age_numeric = parse_numeric_field(age),
  renal_activity_numeric = parse_numeric_field(renal_activity),
  treatment = factor(as.character(treatment)),
  sex = factor(as.character(sex)),
  batch = factor(as.character(batch))
)]
analysis[, timepoint_scaled := as.numeric(scale(timepoint_numeric))]
analysis[, patient_visit_count := .N, by = patient_id]

if (any(!is.finite(analysis$SLEDAI_numeric))) stop("SLEDAI must be numeric for all included SLE visits.")
within_check <- analysis[, .(within_sum = sum(SLEDAI_within)), by = patient_id]
if (max(abs(within_check$within_sum)) > 1e-10) stop("Within-patient SLEDAI centering failed.")

qc_sensitivity_analysis <- data.table::copy(analysis[qc_outlier_flag == FALSE])
center_within_patient(qc_sensitivity_analysis, "SLEDAI", "SLEDAI")
qc_sensitivity_analysis[, timepoint_scaled := as.numeric(scale(timepoint_numeric))]
qc_sensitivity_analysis[, patient_visit_count := .N, by = patient_id]
qc_within_check <- qc_sensitivity_analysis[, .(
  within_sum = sum(SLEDAI_within)
), by = patient_id]
if (max(abs(qc_within_check$within_sum)) > 1e-10) {
  stop("QC sensitivity within-patient SLEDAI centering failed.")
}

programs <- rownames(score_object$score_matrix)
biological_programs <- c("CYTOTOXICITY_SHARED", "CD8_GZMK_STATE")
result_list <- list()
diagnostic_list <- list()
add_model <- function(data, outcome, rhs, family, term) {
  fitted <- fit_target_term(data, outcome, rhs, family, term)
  result_list[[length(result_list) + 1L]] <<- fitted$result
  diagnostic_list[[length(diagnostic_list) + 1L]] <<- fitted$diagnostic
  invisible(fitted)
}

base_rhs <- "SLEDAI_within + SLEDAI_between + timepoint_scaled + (1 | patient_id)"
for (program in programs) {
  add_model(analysis, program, base_rhs, "sle_activity_base", "SLEDAI_within")
  add_model(analysis[patient_visit_count >= 2L], program, base_rhs,
            "sle_activity_repeated_patients", "SLEDAI_within")
}
for (program in programs) {
  add_model(qc_sensitivity_analysis, program, base_rhs,
            "sle_activity_qc_sensitivity", "SLEDAI_within")
}

composition_rhs <- paste(
  "SLEDAI_within + SLEDAI_between + timepoint_scaled +",
  "NK + T_CELL + NEUTROPHIL + ERYTHROID + PLATELET + (1 | patient_id)"
)
clinical_rhs <- paste(
  "SLEDAI_within + SLEDAI_between + timepoint_scaled + renal_activity_numeric +",
  "treatment + age_numeric + sex + batch + (1 | patient_id)"
)
for (program in biological_programs) {
  add_model(analysis, program, composition_rhs, "sle_activity_composition_adjusted", "SLEDAI_within")
  add_model(analysis, program, clinical_rhs, "sle_activity_clinical_adjusted", "SLEDAI_within")
}

secondary_models <- list(
  renal_activity = list(
    rhs = "renal_activity_numeric + SLEDAI_within + SLEDAI_between + timepoint_scaled + (1 | patient_id)",
    term = "renal_activity_numeric"
  ),
  c3 = list(
    rhs = "C3_within + C3_between + timepoint_scaled + (1 | patient_id)",
    term = "C3_within"
  ),
  c4 = list(
    rhs = "C4_within + C4_between + timepoint_scaled + (1 | patient_id)",
    term = "C4_within"
  )
)
for (family in names(secondary_models)) {
  definition <- secondary_models[[family]]
  for (program in programs) {
    add_model(analysis, program, definition$rhs, paste0(family, "_association"), definition$term)
  }
}

model_results <- data.table::rbindlist(result_list, use.names = TRUE, fill = TRUE)
model_results[, fdr := stats::p.adjust(p_value, method = "BH"), by = .(model_family, term)]
data.table::setcolorder(model_results, c(
  "model_family", "program", "term", "estimate", "std_error", "conf_low",
  "conf_high", "p_value", "fdr", "n_patients", "n_visits", "formula",
  "converged", "singular", "warning"
))
diagnostics <- data.table::rbindlist(diagnostic_list, use.names = TRUE, fill = TRUE)

log_message("Running leave-one-patient-out analysis for the primary program.")
primary_program <- cfg$longitudinal$primary_program
patients <- sort(unique(analysis$patient_id))
loo <- data.table::rbindlist(lapply(patients, function(excluded_patient) {
  fitted <- fit_target_term(
    analysis[patient_id != excluded_patient], primary_program, base_rhs,
    "primary_leave_one_patient_out", "SLEDAI_within"
  )$result
  fitted[, excluded_patient := excluded_patient]
  fitted
}), use.names = TRUE, fill = TRUE)
data.table::setcolorder(loo, c("excluded_patient", setdiff(names(loo), "excluded_patient")))

dsdna_values <- trimws(as.character(analysis$anti_dsDNA))
dsdna_format <- data.table::fcase(
  dsdna_values %in% c("", "not_available", "not_applicable", "Data Not Available", "unknown"), "missing",
  grepl("^[0-9]+:[0-9]+$", dsdna_values), "titer",
  is.finite(suppressWarnings(as.numeric(dsdna_values))), "numeric",
  dsdna_values %in% c("(-)", "negative", "Negative"), "negative_text",
  default = "other_text"
)
anti_dsdna_audit <- data.table::data.table(value = dsdna_values, format = dsdna_format)[, .(
  n_visits = .N,
  example_values = paste(head(sort(unique(value)), 8L), collapse = ";")
), by = format]

safe_write_tsv(analysis, outputs[[1]])
safe_write_tsv(model_results, outputs[[2]])
safe_write_tsv(loo, outputs[[3]])
safe_write_tsv(diagnostics, outputs[[4]])
safe_write_tsv(anti_dsdna_audit, outputs[[5]])
atomic_save_rds(list(
  analysis_data = as.data.frame(analysis),
  qc_sensitivity_data = as.data.frame(qc_sensitivity_analysis),
  model_results = as.data.frame(model_results),
  diagnostics = as.data.frame(diagnostics),
  leave_one_patient_out = as.data.frame(loo)
), outputs[[6]], compress = TRUE)
write_session_info(file.path(cfg$paths$logs, "03_longitudinal_models_sessionInfo.txt"))
mark_step_complete(output_dir, signature, outputs)

primary_result <- model_results[
  model_family == "sle_activity_base" & program == primary_program & term == "SLEDAI_within"
]
log_message("Longitudinal models complete. Primary estimate = ",
            signif(primary_result$estimate, 4), ", p = ", signif(primary_result$p_value, 4),
            ", FDR = ", signif(primary_result$fdr, 4), ".")
