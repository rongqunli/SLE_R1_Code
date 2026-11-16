if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")




options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages(library(data.table))
project_root <- "F:/program/Hei/SLE_v1"
output_root <- "F:/program/Hei/SLE_R1/checks/paperwide_fdr_sensitivity_20261115/unassessed14_recheck/heldout_inputs"
source_results <- file.path(project_root, "results/reviewer_response/02_clinical_landscape_validation")
source_script <- file.path(project_root, "analysis/reviewer_clinical_landscape_validation.R")
grouping_path <- file.path(project_root, "pipeline/3.4/data/processed/frozen_grouping_model.rds")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
source_lines <- readLines(source_script, warn = FALSE)
writeLines(sprintf("%04d | %s", seq_along(source_lines), source_lines),
           file.path(output_root, "source_script_numbered.txt"))


grouping_model <- readRDS(grouping_path)
assignments <- as.data.table(grouping_model$assignments)
assignments[, frozen_row := .I]
programs <- grouping_model$features
coordinate_terms <- paste0("coord_", programs)
low_group <- "S02_IMMUNE_LOW"
group_levels <- grouping_model$centroids$subtype
raw_matrix <- as.matrix(assignments[, ..programs])
coordinate_matrix <- sweep(sweep(raw_matrix, 2, grouping_model$center[programs], "-"),
                           2, grouping_model$scale[programs], "/")
assignments[, (coordinate_terms) := as.data.table(coordinate_matrix)]
low_centroid <- as.numeric(grouping_model$centroids[subtype == low_group, ..programs])
distance_raw <- sqrt(rowSums((coordinate_matrix - matrix(low_centroid, nrow(coordinate_matrix),
                                                         length(programs), byrow = TRUE))^2))
distance_center <- mean(distance_raw[assignments$cohort_partition == "discovery"])
distance_scale <- sd(distance_raw[assignments$cohort_partition == "discovery"])
assignments[, landscape_distance_z := (distance_raw - distance_center) / distance_scale]
for (variable in c("SLEDAI", "C3", "C4", "anti_dsDNA", "age")) {
  assignments[, (variable) := suppressWarnings(as.numeric(get(variable)))]
}
assignments[, anti_dsDNA_log := log1p(pmax(anti_dsDNA, 0))]
assignments[, grouping := factor(subtype, levels = group_levels)]
assignments[, high_any := factor(fifelse(subtype == low_group, "IMMUNE_LOW", "PROGRAM_HIGH"),
                                 levels = c("IMMUNE_LOW", "PROGRAM_HIGH"))]
assignments[, race_model := as.character(race)]
assignments[race_model %chin% c("multiple", "native_hw_pac_island"), race_model := "other_or_multiple"]
for (variable in c("sex", "race_model", "batch")) {
  discovery_levels <- sort(unique(assignments[cohort_partition == "discovery" & !is.na(get(variable)), get(variable)]))
  assignments[, (variable) := factor(get(variable), levels = discovery_levels)]
}
clinical_specs <- data.table(trait = c("SLEDAI", "C3", "C4", "anti_dsDNA"),
                             outcome = c("SLEDAI", "C3", "C4", "anti_dsDNA_log"))
covariates <- c("age", "sex", "race_model", "batch")
model_specs <- list(four_group = "grouping", continuous_three_program = coordinate_terms)
make_formula <- function(extra_terms) {
  as.formula(paste("outcome_z ~", paste(c(covariates, extra_terms), collapse = " + ")))
}
prediction_r2 <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]
  if (length(observed) < 5L || sum((observed - mean(observed))^2) <= 0) return(NA_real_)
  1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
}
prediction_rmse <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  sqrt(mean((observed[keep] - predicted[keep])^2))
}
safe_predict <- function(fit, newdata) {
  tryCatch(suppressWarnings(as.numeric(predict(fit, newdata = newdata))),
           error = function(e) rep(NA_real_, nrow(newdata)))
}
prepare_trait <- function(trait, outcome) {
  dat <- copy(assignments)
  dat[, outcome_model := get(outcome)]
  required <- c("outcome_model", covariates, "grouping", "high_any", "landscape_distance_z", coordinate_terms)
  dat <- dat[complete.cases(dat[, ..required]) & is.finite(outcome_model)]
  discovery <- dat[cohort_partition == "discovery"]
  outcome_center <- mean(discovery$outcome_model)
  outcome_scale <- sd(discovery$outcome_model)
  dat[, outcome_z := (outcome_model - outcome_center) / outcome_scale]
  list(all = dat, discovery = dat[cohort_partition == "discovery"],
       outcome_center = outcome_center, outcome_scale = outcome_scale, trait = trait)
}

predictions <- training_membership <- prepared_rows <- normalization <- coefficients <- list()
training_metrics <- heldout_metrics <- comparisons <- ci_checks <- bootstrap_replays <- audits <- list()
recovered_models <- list()
tol <- 1e-12
original_training <- fread(file.path(source_results, "apparent_model_metrics.tsv"))
original_heldout <- fread(file.path(source_results, "trial_linked_heldout_validation_metrics.tsv"))
original_comparisons <- fread(file.path(source_results, "heldout_continuous_vs_four_group.tsv"))
original_flow <- fread(file.path(source_results, "complete_case_flow.tsv"))
stopifnot(nrow(assignments) == 1760L, sum(assignments$cohort_partition == "discovery") == 1319L)
stopifnot(!anyDuplicated(assignments$subject_id), !anyDuplicated(assignments$sample_id))
stopifnot(identical(assignments$subject_id, assignments$patient_id))

for (spec_index in seq_len(nrow(clinical_specs))) {
  trait <- clinical_specs$trait[[spec_index]]
  trait_value <- trait
  outcome <- clinical_specs$outcome[[spec_index]]
  prepared <- prepare_trait(trait, outcome)
  discovery <- prepared$discovery
  all_trait <- prepared$all
  
  full_fits <- lapply(model_specs, function(extra_terms) lm(make_formula(extra_terms), data = discovery))
  recovered_models[[trait]] <- full_fits
  normalization[[trait]] <- data.table(trait, outcome_column = outcome,
                                       transformation = if (trait == "anti_dsDNA") "log1p(pmax(raw, 0))" else "identity",
                                       discovery_n = nrow(discovery), center = prepared$outcome_center,
                                       scale = prepared$outcome_scale)
  training_membership[[trait]] <- discovery[, .(trait = trait, frozen_row, sample_id, subject_id,
                                                trial = cohort_partition, source_trial = as.character(batch))]
  prepared_rows[[trait]] <- copy(all_trait)[, trait := trait]
  for (model_name in names(full_fits)) {
    fit <- full_fits[[model_name]]
    fit_summary <- summary(fit)
    coefficients[[paste(trait, model_name)]] <- data.table(trait, model = model_name,
                                                          term = names(coef(fit)), estimate = unname(coef(fit)))
    old <- original_training[trait == trait_value & model == model_name]
    stopifnot(nrow(old) == 1L)
    training_metrics[[paste(trait, model_name)]] <- data.table(
      trait, model = model_name, n = nobs(fit), old_n = old$n,
      r_squared = fit_summary$r.squared, old_r_squared = old$r_squared,
      r_squared_abs_error = abs(fit_summary$r.squared - old$r_squared),
      adjusted_r_squared = fit_summary$adj.r.squared, old_adjusted_r_squared = old$adjusted_r_squared,
      adjusted_r_squared_abs_error = abs(fit_summary$adj.r.squared - old$adjusted_r_squared),
      parameters = length(coef(fit)), old_parameters = old$parameters,
      formula = paste(deparse(formula(fit)), collapse = " "))
  }
  validation_sets <- list(
    ILLUMINATE1_validation = all_trait[cohort_partition == "ILLUMINATE1_validation"],
    ILLUMINATE2_validation = all_trait[cohort_partition == "ILLUMINATE2_validation"],
    pooled_trial_linked_validation = all_trait[cohort_partition != "discovery"])
  pooled <- validation_sets$pooled_trial_linked_validation
  i1 <- validation_sets$ILLUMINATE1_validation
  i2 <- validation_sets$ILLUMINATE2_validation
  audits[[trait]] <- data.table(
    trait, total_frozen_discovery = sum(assignments$cohort_partition == "discovery"),
    discovery_complete_case_n = nrow(discovery), i1_complete_case_n = nrow(i1), i2_complete_case_n = nrow(i2),
    pooled_complete_case_n = nrow(pooled),
    pooled_is_exact_union = setequal(pooled$subject_id, c(i1$subject_id, i2$subject_id)) && nrow(pooled) == nrow(i1) + nrow(i2),
    trial_overlap_n = length(intersect(i1$subject_id, i2$subject_id)),
    discovery_validation_overlap_n = length(intersect(discovery$subject_id, pooled$subject_id)),
    duplicated_subjects_in_pooled = sum(duplicated(pooled$subject_id)),
    duplicated_samples_in_pooled = sum(duplicated(pooled$sample_id)),
    all_trial_labels_correct = all(as.character(i1$batch) == "ILLUMINATE-1") && all(as.character(i2$batch) == "ILLUMINATE-2"))
  predictions[[trait]] <- data.table(
    trait, frozen_row = pooled$frozen_row, sample_id = pooled$sample_id, subject_id = pooled$subject_id,
    trial = pooled$cohort_partition, source_trial = as.character(pooled$batch),
    observed = pooled$outcome_z,
    pred_continuous = safe_predict(full_fits$continuous_three_program, pooled),
    pred_four_group = safe_predict(full_fits$four_group, pooled),
    outcome_raw = pooled[[trait]], outcome_model = pooled$outcome_model,
    outcome_center = prepared$outcome_center, outcome_scale = prepared$outcome_scale)
  for (validation_name in names(validation_sets)) {
    validation <- validation_sets[[validation_name]]
    observed <- validation$outcome_z
    pred_four <- safe_predict(full_fits$four_group, validation)
    pred_cont <- safe_predict(full_fits$continuous_three_program, validation)
    stopifnot(all(is.finite(observed)), all(is.finite(pred_four)), all(is.finite(pred_cont)))
    r2_four <- prediction_r2(observed, pred_four)
    r2_cont <- prediction_r2(observed, pred_cont)
    for (model_name in names(full_fits)) {
      pred <- if (model_name == "four_group") pred_four else pred_cont
      old <- original_heldout[trait == trait_value & validation_set == validation_name & model == model_name]
      stopifnot(nrow(old) == 1L)
      heldout_metrics[[paste(trait, validation_name, model_name)]] <- data.table(
        trait, validation_set = validation_name, model = model_name, n = length(observed), old_n = old$n,
        r_squared = prediction_r2(observed, pred), old_r_squared = old$heldout_r_squared,
        r_squared_abs_error = abs(prediction_r2(observed, pred) - old$heldout_r_squared),
        rmse = prediction_rmse(observed, pred), old_rmse = old$heldout_rmse,
        rmse_abs_error = abs(prediction_rmse(observed, pred) - old$heldout_rmse))
    }
    old_comp <- original_comparisons[trait == trait_value & validation_set == validation_name]
    stopifnot(nrow(old_comp) == 1L)
    comparisons[[paste(trait, validation_name)]] <- data.table(
      trait, validation_set = validation_name, n = length(observed),
      observed_delta = r2_cont - r2_four,
      old_observed_delta = old_comp$observed_delta_r2_continuous_minus_four,
      delta_abs_error = abs(r2_cont - r2_four - old_comp$observed_delta_r2_continuous_minus_four))
    
    
    seed <- 202608122L + spec_index * 100L + match(validation_name, names(validation_sets))
    set.seed(seed)
    replay <- rbindlist(lapply(seq_len(2000L), function(b) {
      sampled <- sample.int(length(observed), length(observed), replace = TRUE)
      r2f <- prediction_r2(observed[sampled], pred_four[sampled])
      r2c <- prediction_r2(observed[sampled], pred_cont[sampled])
      data.table(bootstrap = b, r2_four_group = r2f, r2_continuous = r2c, delta_continuous_minus_four = r2c - r2f)
    }))
    replay[, `:=`(trait = trait, validation_set = validation_name, seed = seed)]
    bootstrap_replays[[paste(trait, validation_name)]] <- replay
    q <- quantile(replay$delta_continuous_minus_four, c(.025, .975), na.rm = TRUE)
    ci_checks[[paste(trait, validation_name)]] <- data.table(
      trait, validation_set = validation_name, seed, bootstrap_draws = nrow(replay),
      delta_q025 = q[[1]], old_delta_q025 = old_comp$bootstrap_delta_q025,
      q025_abs_error = abs(q[[1]] - old_comp$bootstrap_delta_q025),
      delta_q975 = q[[2]], old_delta_q975 = old_comp$bootstrap_delta_q975,
      q975_abs_error = abs(q[[2]] - old_comp$bootstrap_delta_q975))
  }
}

predictions <- rbindlist(predictions)
training_metrics <- rbindlist(training_metrics)
heldout_metrics <- rbindlist(heldout_metrics)
comparisons <- rbindlist(comparisons)
ci_checks <- rbindlist(ci_checks)
audits <- rbindlist(audits)
all_prepared <- rbindlist(prepared_rows)
flow <- all_prepared[, .(available_complete_cases = .N), by = .(cohort_partition, trait)]
flow_check <- merge(flow, original_flow, by = c("cohort_partition", "trait"), suffixes = c("_restored", "_original"), all = TRUE)
flow_check[, exact_match := available_complete_cases_restored == available_complete_cases_original]
outputs <- list(
  patient_predictions = predictions,
  discovery_training_membership = rbindlist(training_membership),
  normalization = rbindlist(normalization),
  recovered_coefficients = rbindlist(coefficients),
  training_metrics_reproduction = training_metrics,
  heldout_metrics_reproduction = heldout_metrics,
  heldout_delta_reproduction = comparisons,
  legacy_ci_reproduction = ci_checks,
  participant_partition_audit = audits,
  complete_case_reproduction = flow_check)
for (name in names(outputs)) fwrite(outputs[[name]], file.path(output_root, paste0(name, ".tsv")), sep = "\t", na = "NA")
fwrite(rbindlist(bootstrap_replays), file.path(output_root, "legacy_heldout_bootstrap_replay.tsv.gz"), sep = "\t")
saveRDS(list(models = recovered_models, normalization = rbindlist(normalization),
             program_center = grouping_model$center, program_scale = grouping_model$scale,
             group_levels = group_levels, model_specs = model_specs,
             factor_levels = lapply(assignments[, c("sex", "race_model", "batch", "grouping"), with = FALSE], levels),
             provenance = "Restorative deterministic OLS refit on original discovery complete cases; not recovered cache"),
        file.path(output_root, "restored_discovery_clinical_models.rds"))
saveRDS(all_prepared, file.path(output_root, "prepared_trait_complete_cases.rds"))
saveRDS(predictions, file.path(output_root, "patient_predictions.rds"))
capture.output(sessionInfo(), file = file.path(output_root, "sessionInfo.txt"))
stopifnot(nrow(training_metrics) == 8L, nrow(heldout_metrics) == 24L, nrow(comparisons) == 12L,
          all(training_metrics$n == training_metrics$old_n),
          all(training_metrics$parameters == training_metrics$old_parameters),
          max(training_metrics$r_squared_abs_error) < tol,
          max(training_metrics$adjusted_r_squared_abs_error) < tol,
          all(heldout_metrics$n == heldout_metrics$old_n),
          max(heldout_metrics$r_squared_abs_error) < tol,
          max(heldout_metrics$rmse_abs_error) < tol,
          max(comparisons$delta_abs_error) < tol,
          max(ci_checks$q025_abs_error, ci_checks$q975_abs_error) < tol,
          all(flow_check$exact_match), all(audits$pooled_is_exact_union),
          all(audits$trial_overlap_n == 0L), all(audits$discovery_validation_overlap_n == 0L),
          all(audits$duplicated_subjects_in_pooled == 0L), all(audits$duplicated_samples_in_pooled == 0L),
          all(audits$all_trial_labels_correct), !anyDuplicated(predictions[, .(trait, subject_id)]))
summary <- c(
  "PASS: all restoration and participant-partition checks passed at 1e-12 tolerance.",
  sprintf("Recovered fits: %d; held-out model R2/RMSE rows: %d; deltas and legacy CIs: %d.",
          nrow(training_metrics), nrow(heldout_metrics), nrow(comparisons)),
  sprintf("Patient-trait prediction rows: %d; unique held-out participants: %d.", nrow(predictions), uniqueN(predictions$subject_id)),
  sprintf("Maximum training R2 absolute error: %.17g", max(training_metrics$r_squared_abs_error)),
  sprintf("Maximum held-out R2 absolute error: %.17g", max(heldout_metrics$r_squared_abs_error)),
  sprintf("Maximum delta R2 absolute error: %.17g", max(comparisons$delta_abs_error)),
  sprintf("Maximum old CI endpoint absolute error: %.17g", max(ci_checks$q025_abs_error, ci_checks$q975_abs_error)),
  "No validation fitting/tuning, CV rerun, discovery bootstrap rerun, new p-value, or FDR calculation was performed.")
writeLines(summary, file.path(output_root, "RESTORATION_STATUS.txt"))
cat(paste(summary, collapse = "\n"), "\n")
