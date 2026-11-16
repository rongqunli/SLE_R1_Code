if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- if (length(file_arg)) {
  normalizePath(file_arg[[1]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "analysis", "exploratory_program_landscape.R"), winslash = "/", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "exploratory_program_landscape")
for (directory in c("00_audit", "01_pairwise_clinical", "02_continuous_landscape", "03_pharmacodynamics", "figures", "logs")) {
  dir.create(file.path(results_root, directory), recursive = TRUE, showWarnings = FALSE)
}

required_packages <- c("data.table", "ggplot2", "lme4")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

options(stringsAsFactors = FALSE, width = 160)
log_file <- file.path(results_root, "logs", "analysis.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({
  cat("\nSession information\n")
  print(sessionInfo())
  sink()
}, add = TRUE)

primary_programs <- c(
  "CYAN_NEUTROPHIL_IMMATURE",
  "GREENYELLOW_B_CELL_COMPLEMENT",
  "MAGENTA_IFN_ANTIVIRAL"
)
program_labels <- c(
  CYAN_NEUTROPHIL_IMMATURE = "Immature-neutrophil",
  GREENYELLOW_B_CELL_COMPLEMENT = "B-cell/complement",
  MAGENTA_IFN_ANTIVIRAL = "IFN/antiviral",
  landscape_distance = "Distance from IMMUNE_LOW"
)
low_group <- "S02_IMMUNE_LOW"
high_group_to_program <- c(
  S01_NEUTROPHIL_HIGH = "CYAN_NEUTROPHIL_IMMATURE",
  S03_B_CELL_HIGH = "GREENYELLOW_B_CELL_COMPLEMENT",
  S04_IFN_HIGH = "MAGENTA_IFN_ANTIVIRAL"
)
group_labels <- c(
  S01_NEUTROPHIL_HIGH = "NEUTROPHIL_HIGH",
  S02_IMMUNE_LOW = "IMMUNE_LOW",
  S03_B_CELL_HIGH = "B_CELL_HIGH",
  S04_IFN_HIGH = "IFN_HIGH",
  HEALTHY_CONTROL = "HEALTHY_CONTROL"
)

input_paths <- c(
  adult_object = file.path(project_root, "pipeline", "3.3", "data", "processed", "GSE88884_analysis_results.rds"),
  adult_metadata = file.path(project_root, "pipeline", "3.3", "data", "metadata", "GSE88884_meta.tsv"),
  grouping_model = file.path(project_root, "pipeline", "3.4", "data", "processed", "frozen_grouping_model.rds"),
  trial1_scores = file.path(project_root, "pipeline", "3.4", "data", "processed", "ILLUMINATE1_sample_program_scores.tsv"),
  trial2_scores = file.path(project_root, "pipeline", "3.4", "data", "processed", "ILLUMINATE2_sample_program_scores.tsv")
)
if (any(!file.exists(input_paths))) stop("Missing inputs: ", paste(names(input_paths)[!file.exists(input_paths)], collapse = ", "))
input_audit <- data.table::data.table(
  input = names(input_paths),
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  size_bytes = as.numeric(file.info(input_paths)$size),
  md5 = as.character(tools::md5sum(input_paths))
)
data.table::fwrite(input_audit, file.path(results_root, "00_audit", "input_manifest.tsv"), sep = "\t")

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

available_covariates <- function(dat, candidates) {
  candidates[vapply(candidates, function(variable) {
    values <- dat[[variable]]
    values <- values[!is.na(values)]
    length(unique(values)) > 1L
  }, logical(1))]
}

formula_from_terms <- function(outcome, terms) {
  if (!length(terms)) return(stats::as.formula(paste(outcome, "~ 1")))
  stats::as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
}

extract_fixed_term <- function(model, term) {
  coefficients <- coef(summary(model))
  if (!(term %in% rownames(coefficients))) {
    return(data.table::data.table(
      term = term, estimate = NA_real_, std_error = NA_real_, conf_low = NA_real_,
      conf_high = NA_real_, statistic = NA_real_, p_value = NA_real_
    ))
  }
  estimate <- unname(coefficients[term, "Estimate"])
  std_error <- unname(coefficients[term, "Std. Error"])
  statistic <- estimate / std_error
  data.table::data.table(
    term = term,
    estimate = estimate,
    std_error = std_error,
    conf_low = estimate - 1.96 * std_error,
    conf_high = estimate + 1.96 * std_error,
    statistic = statistic,
    p_value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  )
}

fixed_effect_meta <- function(tab, group_columns) {
  tab[is.finite(estimate) & is.finite(std_error) & std_error > 0, {
    weights <- 1 / std_error^2
    estimate_meta <- sum(weights * estimate) / sum(weights)
    standard_error_meta <- sqrt(1 / sum(weights))
    statistic_meta <- estimate_meta / standard_error_meta
    q_statistic <- sum(weights * (estimate - estimate_meta)^2)
    q_df <- .N - 1L
    q_p_value <- if (q_df > 0L) stats::pchisq(q_statistic, df = q_df, lower.tail = FALSE) else NA_real_
    i2 <- if (q_df > 0L && q_statistic > 0) max(0, (q_statistic - q_df) / q_statistic) * 100 else 0
    list(
      trials = .N,
      meta_estimate = estimate_meta,
      meta_std_error = standard_error_meta,
      conf_low = estimate_meta - 1.96 * standard_error_meta,
      conf_high = estimate_meta + 1.96 * standard_error_meta,
      statistic = statistic_meta,
      p_value = 2 * stats::pnorm(abs(statistic_meta), lower.tail = FALSE),
      direction_consistent = length(unique(sign(estimate[estimate != 0]))) <= 1L,
      both_trials_nominal = .N == 2L && all(p_value < 0.05),
      q_statistic = q_statistic,
      q_p_value = q_p_value,
      i2 = i2
    )
  }, by = group_columns]
}

cat("Reading frozen adult and trial inputs\n")
adult_object <- readRDS(input_paths[["adult_object"]])
grouping_model <- readRDS(input_paths[["grouping_model"]])
adult_assignments <- data.table::as.data.table(grouping_model$assignments)
adult_metadata <- data.table::fread(
  input_paths[["adult_metadata"]], sep = "\t",
  na.strings = c("", "NA", "not_available", "not_applicable", "--")
)

if (!all(primary_programs %in% names(adult_assignments))) stop("Adult grouping assignments lack frozen program scores")
if (!all(primary_programs %in% names(grouping_model$centroids))) stop("Grouping model centroids lack frozen program coordinates")

center_values <- grouping_model$center[primary_programs]
scale_values <- grouping_model$scale[primary_programs]
centroid_table <- data.table::as.data.table(grouping_model$centroids)
low_centroid <- as.numeric(centroid_table[subtype == low_group, ..primary_programs])
if (length(low_centroid) != length(primary_programs)) stop("IMMUNE_LOW centroid was not found")

add_landscape_coordinates <- function(tab) {
  result <- data.table::copy(tab)
  raw_matrix <- as.matrix(result[, ..primary_programs])
  coordinate_matrix <- sweep(sweep(raw_matrix, 2, center_values, "-"), 2, scale_values, "/")
  coordinate_names <- paste0("coord_", primary_programs)
  result[, (coordinate_names) := data.table::as.data.table(coordinate_matrix)]
  centered_from_low <- coordinate_matrix - matrix(low_centroid, nrow(coordinate_matrix), ncol(coordinate_matrix), byrow = TRUE)
  result[, landscape_distance := sqrt(rowSums(centered_from_low^2))]
  result[, landscape_distance_z := as.numeric(scale(landscape_distance))]
  result
}

adult_sle <- add_landscape_coordinates(adult_assignments)
adult_sle[, grouping := factor(subtype, levels = names(group_labels)[1:4])]
adult_sle[, high_any := factor(ifelse(subtype == low_group, "IMMUNE_LOW", "PROGRAM_HIGH"), levels = c("IMMUNE_LOW", "PROGRAM_HIGH"))]
for (variable in c("SLEDAI", "C3", "C4", "anti_dsDNA", "age")) adult_sle[, (variable) := safe_numeric(get(variable))]
adult_sle[, anti_dsDNA_log := log1p(pmax(anti_dsDNA, 0))]
for (variable in c("sex", "race", "batch")) adult_sle[, (variable) := factor(get(variable))]

discovery_sle <- adult_sle[cohort_partition == "discovery"]
if (nrow(discovery_sle) != 1319L) warning("Expected 1,319 discovery SLE subjects; observed ", nrow(discovery_sle))

clinical_specs <- data.table::data.table(
  trait = c("SLEDAI", "C3", "C4", "anti_dsDNA", "age"),
  model_variable = c("SLEDAI", "C3", "C4", "anti_dsDNA_log", "age"),
  activity_direction = c(1, -1, -1, 1, NA_real_)
)

clinical_group_summary <- data.table::rbindlist(lapply(clinical_specs$trait, function(trait_name) {
  discovery_sle[, .(
    n = sum(is.finite(get(trait_name))),
    median = stats::median(get(trait_name), na.rm = TRUE),
    q25 = stats::quantile(get(trait_name), 0.25, na.rm = TRUE),
    q75 = stats::quantile(get(trait_name), 0.75, na.rm = TRUE)
  ), by = subtype][, trait := trait_name]
}), fill = TRUE)
data.table::fwrite(clinical_group_summary, file.path(results_root, "01_pairwise_clinical", "clinical_group_summary.tsv"), sep = "\t")

fit_pairwise_clinical <- function(high_group, trait_name, model_variable) {
  dat <- data.table::copy(discovery_sle[subtype %in% c(low_group, high_group)])
  dat[, high_state := as.integer(subtype == high_group)]
  dat[, outcome_raw := get(trait_name)]
  dat[, outcome_model := get(model_variable)]
  candidate_covariates <- if (trait_name == "age") c("sex", "race", "batch") else c("age", "sex", "race", "batch")
  covariates <- available_covariates(dat, candidate_covariates)
  required_columns <- c("outcome_raw", "outcome_model", "high_state", covariates)
  dat <- dat[stats::complete.cases(dat[, ..required_columns]) & is.finite(outcome_model)]
  if (nrow(dat) < 30L || length(unique(dat$high_state)) < 2L) return(NULL)
  dat[, outcome_z := as.numeric(scale(outcome_model))]
  fit <- stats::lm(formula_from_terms("outcome_z", c("high_state", covariates)), data = dat)
  coefficient_table <- coef(summary(fit))
  if (!("high_state" %in% rownames(coefficient_table))) return(NULL)
  estimate <- unname(coefficient_table["high_state", "Estimate"])
  standard_error <- unname(coefficient_table["high_state", "Std. Error"])
  statistic <- unname(coefficient_table["high_state", "t value"])
  p_value <- unname(coefficient_table["high_state", "Pr(>|t|)"])
  critical <- stats::qt(0.975, df = stats::df.residual(fit))
  data.table::data.table(
    high_group = high_group,
    reference_group = low_group,
    trait = trait_name,
    model_variable = model_variable,
    adjusted_standardized_effect = estimate,
    std_error = standard_error,
    conf_low = estimate - critical * standard_error,
    conf_high = estimate + critical * standard_error,
    statistic = statistic,
    p_value = p_value,
    n_high = sum(dat$high_state == 1L),
    n_low = sum(dat$high_state == 0L),
    median_high = stats::median(dat[high_state == 1L, outcome_raw], na.rm = TRUE),
    median_low = stats::median(dat[high_state == 0L, outcome_raw], na.rm = TRUE),
    covariates = paste(covariates, collapse = ";")
  )
}

pairwise_clinical <- data.table::rbindlist(lapply(names(high_group_to_program), function(high_group) {
  data.table::rbindlist(lapply(seq_len(nrow(clinical_specs)), function(index) {
    fit_pairwise_clinical(high_group, clinical_specs$trait[[index]], clinical_specs$model_variable[[index]])
  }), fill = TRUE)
}), fill = TRUE)
pairwise_clinical[, fdr := stats::p.adjust(p_value, method = "BH")]
pairwise_clinical <- merge(pairwise_clinical, clinical_specs[, .(trait, activity_direction)], by = "trait", all.x = TRUE)
pairwise_clinical[, activity_aligned_effect := adjusted_standardized_effect * activity_direction]
data.table::fwrite(pairwise_clinical, file.path(results_root, "01_pairwise_clinical", "program_high_vs_immune_low_clinical_contrasts.tsv"), sep = "\t")

clinical_separation_summary <- pairwise_clinical[trait != "age", .(
  disease_traits = .N,
  fdr_significant_traits = sum(fdr < 0.05, na.rm = TRUE),
  mean_absolute_standardized_effect = mean(abs(adjusted_standardized_effect), na.rm = TRUE),
  mean_activity_aligned_effect = mean(activity_aligned_effect, na.rm = TRUE),
  minimum_fdr = min(fdr, na.rm = TRUE)
), by = high_group]
data.table::setorder(clinical_separation_summary, -fdr_significant_traits, -mean_absolute_standardized_effect)
clinical_separation_summary[, exploratory_rank := seq_len(.N)]
data.table::fwrite(clinical_separation_summary, file.path(results_root, "01_pairwise_clinical", "clinical_separation_ranking.tsv"), sep = "\t")

auc_from_ranks <- function(outcome, score) {
  keep <- is.finite(outcome) & is.finite(score)
  outcome <- outcome[keep]
  score <- score[keep]
  positives <- sum(outcome == 1L)
  negatives <- sum(outcome == 0L)
  if (!positives || !negatives) return(NA_real_)
  score_ranks <- rank(score, ties.method = "average")
  (sum(score_ranks[outcome == 1L]) - positives * (positives + 1) / 2) / (positives * negatives)
}

cross_validated_auc <- function(dat, predictor_terms, repeats = 25L, folds = 10L, seed = 20260807L) {
  required_columns <- c("outcome", predictor_terms)
  work <- as.data.frame(dat[stats::complete.cases(dat[, ..required_columns]), ..required_columns])
  auc_values <- rep(NA_real_, repeats)
  for (repeat_index in seq_len(repeats)) {
    set.seed(seed + repeat_index)
    fold_id <- integer(nrow(work))
    for (class_value in c(0L, 1L)) {
      class_rows <- which(work$outcome == class_value)
      fold_id[class_rows] <- sample(rep(seq_len(folds), length.out = length(class_rows)))
    }
    predictions <- rep(NA_real_, nrow(work))
    for (fold_index in seq_len(folds)) {
      training <- work[fold_id != fold_index, , drop = FALSE]
      testing <- work[fold_id == fold_index, , drop = FALSE]
      fit <- suppressWarnings(tryCatch(
        stats::glm(
          formula_from_terms("outcome", predictor_terms),
          data = training,
          family = stats::binomial(),
          control = stats::glm.control(maxit = 100)
        ),
        error = function(error) NULL
      ))
      if (!is.null(fit) && nrow(testing)) {
        predictions[fold_id == fold_index] <- suppressWarnings(stats::predict(fit, newdata = testing, type = "response"))
      }
    }
    auc_values[[repeat_index]] <- auc_from_ranks(work$outcome, predictions)
  }
  auc_values <- auc_values[is.finite(auc_values)]
  data.table::data.table(
    n = nrow(work),
    positives = sum(work$outcome == 1L),
    negatives = sum(work$outcome == 0L),
    repeats_completed = length(auc_values),
    folds = folds,
    mean_auc = mean(auc_values),
    sd_auc = stats::sd(auc_values),
    empirical_low = stats::quantile(auc_values, 0.025, names = FALSE),
    empirical_high = stats::quantile(auc_values, 0.975, names = FALSE)
  )
}

coordinate_terms <- paste0("coord_", primary_programs)
pairwise_molecular_discrimination <- data.table::rbindlist(Map(function(high_group, seed_offset) {
  dat <- data.table::copy(discovery_sle[subtype %in% c(low_group, high_group)])
  dat[, outcome := as.integer(subtype == high_group)]
  focal_term <- paste0("coord_", high_group_to_program[[high_group]])
  model_specs <- list(
    focal_program = focal_term,
    distance_from_immune_low = "landscape_distance_z",
    continuous_three_program = coordinate_terms
  )
  data.table::rbindlist(Map(function(terms, model_name, model_offset) {
    cross_validated_auc(dat, terms, seed = 20260807L + seed_offset * 1000L + model_offset)[, `:=`(
      comparison = paste(high_group, "vs", low_group),
      high_group = high_group,
      reference_group = low_group,
      model = model_name,
      predictors = paste(terms, collapse = ";")
    )]
  }, model_specs, names(model_specs), seq_along(model_specs)), fill = TRUE)
}, names(high_group_to_program), seq_along(high_group_to_program)), fill = TRUE)

pooled_discrimination_data <- data.table::copy(discovery_sle)
pooled_discrimination_data[, outcome := as.integer(subtype != low_group)]
pooled_model_specs <- list(
  distance_from_immune_low = "landscape_distance_z",
  continuous_three_program = coordinate_terms
)
pooled_molecular_discrimination <- data.table::rbindlist(Map(function(terms, model_name, model_offset) {
  cross_validated_auc(pooled_discrimination_data, terms, seed = 20263807L + model_offset)[, `:=`(
    comparison = paste("ANY_PROGRAM_HIGH", "vs", low_group),
    high_group = "ANY_PROGRAM_HIGH",
    reference_group = low_group,
    model = model_name,
    predictors = paste(terms, collapse = ";")
  )]
}, pooled_model_specs, names(pooled_model_specs), seq_along(pooled_model_specs)), fill = TRUE)

molecular_discrimination <- data.table::rbindlist(
  list(pairwise_molecular_discrimination, pooled_molecular_discrimination),
  use.names = TRUE,
  fill = TRUE
)
molecular_discrimination[, interpretation_note := "Labels were defined from the same frozen program coordinates; cross-validated AUC measures boundary reproducibility, not independent biological validation."]
data.table::setcolorder(
  molecular_discrimination,
  c("comparison", "high_group", "reference_group", "model", "predictors", "n", "positives", "negatives", "repeats_completed", "folds", "mean_auc", "sd_auc", "empirical_low", "empirical_high", "interpretation_note")
)
data.table::fwrite(
  molecular_discrimination,
  file.path(results_root, "02_continuous_landscape", "cross_validated_molecular_discrimination.tsv"),
  sep = "\t"
)

model_metrics <- function(fit, base_fit, trait_name, model_name) {
  model_summary <- summary(fit)
  comparison <- tryCatch(stats::anova(base_fit, fit), error = function(error) NULL)
  incremental_p <- if (!is.null(comparison) && nrow(comparison) >= 2L) comparison$`Pr(>F)`[[2]] else NA_real_
  data.table::data.table(
    trait = trait_name,
    model = model_name,
    n = stats::nobs(fit),
    parameters = length(stats::coef(fit)),
    r_squared = model_summary$r.squared,
    adjusted_r_squared = model_summary$adj.r.squared,
    incremental_r_squared = model_summary$r.squared - summary(base_fit)$r.squared,
    aic = stats::AIC(fit),
    bic = stats::BIC(fit),
    incremental_p_value = incremental_p
  )
}

fit_landscape_model_comparison <- function(trait_name, model_variable) {
  dat <- data.table::copy(discovery_sle)
  dat[, outcome_model := get(model_variable)]
  candidate_covariates <- if (trait_name == "age") c("sex", "race", "batch") else c("age", "sex", "race", "batch")
  coordinate_terms <- paste0("coord_", primary_programs)
  required_columns <- c("outcome_model", candidate_covariates, "grouping", "high_any", "landscape_distance_z", coordinate_terms)
  dat <- dat[stats::complete.cases(dat[, ..required_columns]) & is.finite(outcome_model)]
  if (nrow(dat) < 100L) return(NULL)
  dat[, outcome_z := as.numeric(scale(outcome_model))]
  covariates <- available_covariates(dat, candidate_covariates)
  base_fit <- stats::lm(formula_from_terms("outcome_z", covariates), data = dat)
  fits <- list(
    binary_high_vs_low = stats::lm(formula_from_terms("outcome_z", c(covariates, "high_any")), data = dat),
    four_group = stats::lm(formula_from_terms("outcome_z", c(covariates, "grouping")), data = dat),
    distance_from_immune_low = stats::lm(formula_from_terms("outcome_z", c(covariates, "landscape_distance_z")), data = dat),
    continuous_three_program = stats::lm(formula_from_terms("outcome_z", c(covariates, coordinate_terms)), data = dat)
  )
  data.table::rbindlist(c(
    list(model_metrics(base_fit, base_fit, trait_name, "covariates_only")),
    Map(function(fit, name) model_metrics(fit, base_fit, trait_name, name), fits, names(fits))
  ), fill = TRUE)
}

landscape_model_comparison <- data.table::rbindlist(lapply(seq_len(nrow(clinical_specs)), function(index) {
  fit_landscape_model_comparison(clinical_specs$trait[[index]], clinical_specs$model_variable[[index]])
}), fill = TRUE)
landscape_model_comparison[, incremental_fdr := stats::p.adjust(incremental_p_value, method = "BH")]
landscape_model_comparison[, delta_aic_from_best := aic - min(aic), by = trait]
data.table::fwrite(landscape_model_comparison, file.path(results_root, "02_continuous_landscape", "clinical_model_comparison.tsv"), sep = "\t")

best_landscape_models <- landscape_model_comparison[model != "covariates_only"][order(trait, aic), .SD[1], by = trait]
data.table::fwrite(best_landscape_models, file.path(results_root, "02_continuous_landscape", "best_clinical_representation_by_trait.tsv"), sep = "\t")

adult_score_long <- data.table::as.data.table(adult_object$program_scores)[program %in% primary_programs]
adult_score_wide <- data.table::dcast(
  adult_score_long,
  sample_id + subject_id + condition ~ program,
  value.var = "score_z"
)
adult_score_wide <- add_landscape_coordinates(adult_score_wide)
adult_state_map <- unique(adult_assignments[, .(sample_id, subtype)])
adult_score_wide <- merge(adult_score_wide, adult_state_map, by = "sample_id", all.x = TRUE)
adult_score_wide[, state := ifelse(condition == "SLE", subtype, "HEALTHY_CONTROL")]
adult_score_wide <- adult_score_wide[state %in% names(group_labels)]

coordinate_columns <- paste0("coord_", primary_programs)
state_profiles <- adult_score_wide[, c(
  list(n = .N, mean_landscape_distance = mean(landscape_distance), median_landscape_distance = stats::median(landscape_distance)),
  lapply(.SD, mean)
), by = state, .SDcols = coordinate_columns]
data.table::fwrite(state_profiles, file.path(results_root, "02_continuous_landscape", "adult_state_program_profiles.tsv"), sep = "\t")

pca_fit <- stats::prcomp(as.matrix(adult_score_wide[, ..coordinate_columns]), center = FALSE, scale. = FALSE)
pca_coordinates <- data.table::data.table(
  sample_id = adult_score_wide$sample_id,
  subject_id = adult_score_wide$subject_id,
  state = adult_score_wide$state,
  landscape_distance = adult_score_wide$landscape_distance,
  PC1 = pca_fit$x[, 1],
  PC2 = pca_fit$x[, 2]
)
data.table::fwrite(pca_coordinates, file.path(results_root, "02_continuous_landscape", "adult_program_landscape_pca.tsv"), sep = "\t")

sink()
