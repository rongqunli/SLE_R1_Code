if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
parse_numeric <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "not_available", "not_applicable", "Data_not_available")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

read_matrix_sample_ids <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "!Sample_geo_accession")) {
      return(gsub('^"|"$', "", strsplit(line, "\t", fixed = TRUE)[[1]][-1L]))
    }
  }
  stop("Could not find GEO sample IDs in series matrix.")
}

read_platform_mapping <- function(path) {
  annotation <- data.table::fread(path, skip = "ID\tGene title", sep = "\t", quote = "")
  required <- c("ID", "Gene symbol")
  if (!all(required %in% names(annotation))) stop("GPL annotation columns are incomplete.")
  mapping <- annotation[, .(probe_id = as.character(ID), gene_symbol_raw = as.character(`Gene symbol`))]
  mapping[, mapping_status := data.table::fcase(
    !nzchar(gene_symbol_raw) | gene_symbol_raw == "---", "unmapped",
    grepl(" /// ", gene_symbol_raw, fixed = TRUE), "multi_symbol",
    default = "unique_symbol"
  )]
  mapping[, gene_symbol := ifelse(
    mapping_status == "unique_symbol", toupper(trimws(gene_symbol_raw)), NA_character_
  )]
  unique(mapping, by = "probe_id")
}

center_within_patient <- function(data, source, prefix) {
  numeric_col <- paste0(prefix, "_numeric")
  between_col <- paste0(prefix, "_between")
  within_col <- paste0(prefix, "_within")
  data[, (numeric_col) := parse_numeric(get(source))]
  data[, (between_col) := {
    values <- get(numeric_col)
    if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
  }, by = patient_id]
  data[, (within_col) := get(numeric_col) - get(between_col)]
  invisible(data)
}

sanitize_message <- function(x) {
  x <- paste(unique(x[nzchar(x)]), collapse = " | ")
  gsub("[\r\n\t]+", " ", x)
}

fit_mixed_target <- function(data, outcome, rhs, family, target_term, feature, role) {
  model_data <- data.table::copy(data)
  model_data[, outcome_value := as.numeric(get(outcome))]
  formula_text <- paste("outcome_value ~", rhs)
  formula <- stats::as.formula(formula_text)
  variables <- all.vars(formula)
  model_data <- droplevels(as.data.frame(model_data[stats::complete.cases(model_data[, ..variables])]))
  warnings <- character(0)
  error_message <- ""
  fit <- tryCatch(
    withCallingHandlers(
      lme4::lmer(
        formula, data = model_data, REML = FALSE,
        control = lme4::lmerControl(
          optimizer = "bobyqa", optCtrl = list(maxfun = 200000),
          check.conv.singular = "ignore"
        )
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  result <- data.table::data.table(
    model_family = family, analysis_role = role, feature = feature, term = target_term,
    estimate = NA_real_, std_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
    p_value = NA_real_, n_patients = data.table::uniqueN(model_data$patient_id),
    n_visits = nrow(model_data), formula = formula_text, converged = FALSE,
    singular = NA, warning = sanitize_message(c(warnings, error_message))
  )
  if (is.null(fit)) return(result)
  coefficients <- stats::coef(summary(fit))
  convergence <- fit@optinfo$conv$lme4$messages
  result[, `:=`(
    converged = is.null(convergence),
    singular = lme4::isSingular(fit, tol = 1e-4),
    warning = sanitize_message(c(warnings, convergence))
  )]
  if (!target_term %in% rownames(coefficients)) {
    result[, warning := sanitize_message(c(warning, "Target term absent from fitted coefficients."))]
    return(result)
  }
  beta_value <- unname(coefficients[target_term, "Estimate"])
  standard_error <- unname(coefficients[target_term, "Std. Error"])
  result[, `:=`(
    estimate = beta_value,
    std_error = standard_error,
    conf_low = beta_value - 1.96 * standard_error,
    conf_high = beta_value + 1.96 * standard_error,
    p_value = 2 * stats::pnorm(abs(beta_value / standard_error), lower.tail = FALSE)
  )]
  result
}

signed_z <- function(estimate, p_value) {
  if (!is.finite(estimate) || !is.finite(p_value) || p_value <= 0 || p_value > 1) return(NA_real_)
  sign(estimate) * stats::qnorm(p_value / 2, lower.tail = FALSE)
}
