if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")


options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages(library(data.table))

project_root <- "F:/program/Hei/SLE_v1"
output_root <- file.path(
  project_root, "results", "reviewer_response", "02_clinical_landscape_validation"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

grouping_path <- file.path(
  project_root, "pipeline", "3.4", "data", "processed", "frozen_grouping_model.rds"
)
grouping_model <- readRDS(grouping_path)
assignments <- as.data.table(grouping_model$assignments)

programs <- grouping_model$features
coordinate_terms <- paste0("coord_", programs)
low_group <- "S02_IMMUNE_LOW"
group_levels <- grouping_model$centroids$subtype

raw_matrix <- as.matrix(assignments[, ..programs])
coordinate_matrix <- sweep(
  sweep(raw_matrix, 2, grouping_model$center[programs], "-"),
  2, grouping_model$scale[programs], "/"
)
assignments[, (coordinate_terms) := as.data.table(coordinate_matrix)]
low_centroid <- as.numeric(
  grouping_model$centroids[subtype == low_group, ..programs]
)
distance_raw <- sqrt(rowSums(
  (coordinate_matrix - matrix(low_centroid, nrow(coordinate_matrix), length(programs), byrow = TRUE))^2
))
distance_center <- mean(distance_raw[assignments$cohort_partition == "discovery"])
distance_scale <- sd(distance_raw[assignments$cohort_partition == "discovery"])
assignments[, landscape_distance_z := (distance_raw - distance_center) / distance_scale]

for (variable in c("SLEDAI", "C3", "C4", "anti_dsDNA", "age")) {
  assignments[, (variable) := suppressWarnings(as.numeric(get(variable)))]
}
assignments[, anti_dsDNA_log := log1p(pmax(anti_dsDNA, 0))]
assignments[, grouping := factor(subtype, levels = group_levels)]
assignments[, high_any := factor(
  fifelse(subtype == low_group, "IMMUNE_LOW", "PROGRAM_HIGH"),
  levels = c("IMMUNE_LOW", "PROGRAM_HIGH")
)]



assignments[, race_model := as.character(race)]
assignments[race_model %chin% c("multiple", "native_hw_pac_island"), race_model := "other_or_multiple"]
for (variable in c("sex", "race_model", "batch")) {
  discovery_levels <- sort(unique(assignments[cohort_partition == "discovery" & !is.na(get(variable)), get(variable)]))
  assignments[, (variable) := factor(get(variable), levels = discovery_levels)]
}

clinical_specs <- data.table(
  trait = c("SLEDAI", "C3", "C4", "anti_dsDNA"),
  outcome = c("SLEDAI", "C3", "C4", "anti_dsDNA_log")
)
covariates <- c("age", "sex", "race_model", "batch")
model_specs <- list(
  covariates_only = character(),
  binary_high_vs_low = "high_any",
  distance_from_immune_low = "landscape_distance_z",
  four_group = "grouping",
  continuous_three_program = coordinate_terms
)
molecular_df <- c(
  covariates_only = 0L,
  binary_high_vs_low = 1L,
  distance_from_immune_low = 1L,
  four_group = 3L,
  continuous_three_program = 3L
)

model_definitions <- data.table(
  model = names(model_specs),
  molecular_predictors = vapply(model_specs, function(x) paste(x, collapse = ";"), character(1)),
  added_degrees_of_freedom = unname(molecular_df[names(model_specs)]),
  definition = c(
    "age + sex + race + batch",
    "covariates plus PROGRAM_HIGH versus IMMUNE_LOW indicator",
    "covariates plus discovery-standardized Euclidean distance from the frozen IMMUNE_LOW centroid",
    "covariates plus the four-level nearest-centroid grouping",
    "covariates plus all three discovery-standardized frozen program coordinates"
  )
)
fwrite(model_definitions, file.path(output_root, "model_definitions.tsv"), sep = "\t")

grouping_definition <- data.table(
  item = c(
    "input_variables", "standardization", "consensus_algorithm", "distance",
    "resampling", "candidate_K", "K_selection", "selected_K", "reference_regions",
    "boundary_assignment", "low_margin"
  ),
  implementation = c(
    paste(programs, collapse = ";"),
    "discovery-set mean and standard deviation for each frozen program coordinate",
    "subsampled k-means followed by average-linkage clustering of 1-consensus",
    "Euclidean distance for k-means and squared Euclidean distance to frozen centroids",
    "500 resamples; 80% of discovery participants; k-means nstart=30; iter.max=100",
    "K=2 through K=6",
    "minimum joint rank of PAC and negative mean silhouette among solutions with mean silhouette >=0.20",
    as.character(grouping_model$selected_k),
    "Voronoi regions induced by the four frozen centroids in the standardized three-program space",
    "nearest squared-Euclidean centroid; exact ties resolved by frozen centroid order",
    "assignment margin (d2-d1)/d2 <0.10"
  )
)
fwrite(grouping_definition, file.path(output_root, "grouping_definition.tsv"), sep = "\t")

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
  tryCatch(
    suppressWarnings(as.numeric(predict(fit, newdata = newdata))),
    error = function(e) rep(NA_real_, nrow(newdata))
  )
}

fit_all_models <- function(data) {
  lapply(model_specs, function(extra_terms) lm(make_formula(extra_terms), data = data))
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
  list(
    all = dat, discovery = dat[cohort_partition == "discovery"],
    outcome_center = outcome_center, outcome_scale = outcome_scale,
    trait = trait
  )
}

complete_case_rows <- list()
apparent_rows <- list()
bootstrap_summary_rows <- list()
bootstrap_distribution_rows <- list()
cv_rows <- list()
cv_fold_rows <- list()
cv_comparison_rows <- list()
heldout_rows <- list()
heldout_comparison_rows <- list()

n_bootstrap <- 1000L
n_cv_repeats <- 20L
n_folds <- 10L
n_heldout_bootstrap <- 2000L

for (spec_index in seq_len(nrow(clinical_specs))) {
  trait <- clinical_specs$trait[[spec_index]]
  outcome <- clinical_specs$outcome[[spec_index]]
  prepared <- prepare_trait(trait, outcome)
  discovery <- prepared$discovery
  all_trait <- prepared$all

  complete_case_rows[[trait]] <- all_trait[, .(
    available_complete_cases = .N
  ), by = cohort_partition][, trait := trait]

  full_fits <- fit_all_models(discovery)
  apparent <- rbindlist(lapply(names(full_fits), function(model_name) {
    fit <- full_fits[[model_name]]
    fit_summary <- summary(fit)
    data.table(
      trait = trait, model = model_name, n = nobs(fit),
      parameters = length(coef(fit)), added_degrees_of_freedom = molecular_df[[model_name]],
      r_squared = fit_summary$r.squared,
      adjusted_r_squared = fit_summary$adj.r.squared,
      aic = AIC(fit), bic = BIC(fit)
    )
  }))
  base_r2 <- apparent[model == "covariates_only", r_squared]
  apparent[, incremental_r_squared := r_squared - base_r2]
  apparent_rows[[trait]] <- apparent

  
  trait_cv <- list()
  trait_cv_folds <- list()
  for (repeat_index in seq_len(n_cv_repeats)) {
    set.seed(202608120L + spec_index * 1000L + repeat_index)
    fold_id <- sample(rep(seq_len(n_folds), length.out = nrow(discovery)))
    predictions <- setNames(
      lapply(names(model_specs), function(x) rep(NA_real_, nrow(discovery))),
      names(model_specs)
    )
    for (fold in seq_len(n_folds)) {
      train <- discovery[fold_id != fold]
      test <- discovery[fold_id == fold]
      fold_fits <- fit_all_models(train)
      fold_mse <- numeric(length(fold_fits))
      names(fold_mse) <- names(fold_fits)
      for (model_name in names(fold_fits)) {
        pred <- safe_predict(fold_fits[[model_name]], test)
        predictions[[model_name]][fold_id == fold] <- pred
        fold_mse[[model_name]] <- mean((test$outcome_z - pred)^2, na.rm = TRUE)
      }
      trait_cv_folds[[paste(repeat_index, fold, sep = "_")]] <- data.table(
        trait = trait, cv_repeat = repeat_index, fold = fold,
        mse_four_group = fold_mse[["four_group"]],
        mse_continuous = fold_mse[["continuous_three_program"]],
        mse_difference_four_minus_continuous =
          fold_mse[["four_group"]] - fold_mse[["continuous_three_program"]]
      )
    }
    repeat_metrics <- rbindlist(lapply(names(predictions), function(model_name) {
      data.table(
        trait = trait, cv_repeat = repeat_index, model = model_name,
        cv_r_squared = prediction_r2(discovery$outcome_z, predictions[[model_name]]),
        cv_rmse = prediction_rmse(discovery$outcome_z, predictions[[model_name]])
      )
    }))
    repeat_base <- repeat_metrics[model == "covariates_only", cv_r_squared]
    repeat_metrics[, cv_incremental_r_squared := cv_r_squared - repeat_base]
    trait_cv[[repeat_index]] <- repeat_metrics
  }
  trait_cv <- rbindlist(trait_cv)
  trait_cv_folds <- rbindlist(trait_cv_folds)
  cv_rows[[trait]] <- trait_cv
  cv_fold_rows[[trait]] <- trait_cv_folds

  cv_wide <- dcast(trait_cv, trait + cv_repeat ~ model, value.var = "cv_r_squared")
  cv_wide[, delta_continuous_minus_four := continuous_three_program - four_group]
  fold_difference <- trait_cv_folds$mse_difference_four_minus_continuous
  correction <- 1 / nrow(trait_cv_folds) + 1 / (n_folds - 1)
  corrected_se <- sqrt(correction * var(fold_difference))
  corrected_t <- mean(fold_difference) / corrected_se
  cv_comparison_rows[[trait]] <- data.table(
    trait = trait,
    mean_cv_r2_difference_continuous_minus_four = mean(cv_wide$delta_continuous_minus_four),
    cv_difference_q025 = quantile(cv_wide$delta_continuous_minus_four, 0.025),
    cv_difference_q975 = quantile(cv_wide$delta_continuous_minus_four, 0.975),
    corrected_resampled_t = corrected_t,
    corrected_resampled_p = 2 * pt(abs(corrected_t), df = nrow(trait_cv_folds) - 1L, lower.tail = FALSE),
    repeats = n_cv_repeats, folds = n_folds
  )

  
  trait_bootstrap <- vector("list", n_bootstrap)
  set.seed(202608121L + spec_index)
  for (bootstrap_index in seq_len(n_bootstrap)) {
    sampled <- sample.int(nrow(discovery), nrow(discovery), replace = TRUE)
    boot_data <- discovery[sampled]
    boot_fits <- fit_all_models(boot_data)
    boot_metrics <- rbindlist(lapply(names(boot_fits), function(model_name) {
      fit <- boot_fits[[model_name]]
      data.table(
        bootstrap = bootstrap_index, model = model_name,
        bootstrap_apparent_r2 = prediction_r2(boot_data$outcome_z, safe_predict(fit, boot_data)),
        original_test_r2 = prediction_r2(discovery$outcome_z, safe_predict(fit, discovery))
      )
    }))
    base_app <- boot_metrics[model == "covariates_only", bootstrap_apparent_r2]
    base_test <- boot_metrics[model == "covariates_only", original_test_r2]
    boot_metrics[, `:=`(
      bootstrap_apparent_incremental_r2 = bootstrap_apparent_r2 - base_app,
      original_test_incremental_r2 = original_test_r2 - base_test
    )]
    trait_bootstrap[[bootstrap_index]] <- boot_metrics
  }
  trait_bootstrap <- rbindlist(trait_bootstrap)
  trait_bootstrap[, trait := trait]
  bootstrap_distribution_rows[[trait]] <- trait_bootstrap

  trait_bootstrap_summary <- trait_bootstrap[, .(
    r2_bootstrap_q025 = quantile(bootstrap_apparent_r2, 0.025, na.rm = TRUE),
    r2_bootstrap_q975 = quantile(bootstrap_apparent_r2, 0.975, na.rm = TRUE),
    incremental_r2_bootstrap_q025 = quantile(bootstrap_apparent_incremental_r2, 0.025, na.rm = TRUE),
    incremental_r2_bootstrap_q975 = quantile(bootstrap_apparent_incremental_r2, 0.975, na.rm = TRUE),
    mean_r2_optimism = mean(bootstrap_apparent_r2 - original_test_r2, na.rm = TRUE),
    mean_incremental_r2_optimism = mean(
      bootstrap_apparent_incremental_r2 - original_test_incremental_r2, na.rm = TRUE
    )
  ), by = model]
  trait_bootstrap_summary[, trait := trait]
  trait_bootstrap_summary <- merge(
    trait_bootstrap_summary,
    apparent[, .(trait, model, r_squared, adjusted_r_squared, incremental_r_squared)],
    by = c("trait", "model")
  )
  trait_bootstrap_summary[, `:=`(
    optimism_corrected_r_squared = r_squared - mean_r2_optimism,
    optimism_corrected_incremental_r_squared =
      incremental_r_squared - mean_incremental_r2_optimism
  )]
  bootstrap_summary_rows[[trait]] <- trait_bootstrap_summary

  
  validation_sets <- list(
    ILLUMINATE1_validation = all_trait[cohort_partition == "ILLUMINATE1_validation"],
    ILLUMINATE2_validation = all_trait[cohort_partition == "ILLUMINATE2_validation"],
    pooled_trial_linked_validation = all_trait[cohort_partition != "discovery"]
  )
  for (validation_name in names(validation_sets)) {
    validation <- validation_sets[[validation_name]]
    prediction_table <- data.table(observed = validation$outcome_z)
    for (model_name in names(full_fits)) {
      prediction_table[, (model_name) := safe_predict(full_fits[[model_name]], validation)]
    }
    validation_metrics <- rbindlist(lapply(names(full_fits), function(model_name) {
      data.table(
        trait = trait, validation_set = validation_name, model = model_name,
        n = sum(is.finite(prediction_table$observed) & is.finite(prediction_table[[model_name]])),
        heldout_r_squared = prediction_r2(prediction_table$observed, prediction_table[[model_name]]),
        heldout_rmse = prediction_rmse(prediction_table$observed, prediction_table[[model_name]])
      )
    }))
    heldout_base <- validation_metrics[model == "covariates_only", heldout_r_squared]
    validation_metrics[, heldout_incremental_r_squared := heldout_r_squared - heldout_base]

    set.seed(202608122L + spec_index * 100L + match(validation_name, names(validation_sets)))
    heldout_bootstrap <- vector("list", n_heldout_bootstrap)
    for (bootstrap_index in seq_len(n_heldout_bootstrap)) {
      sampled <- sample.int(nrow(prediction_table), nrow(prediction_table), replace = TRUE)
      sampled_predictions <- prediction_table[sampled]
      model_r2 <- vapply(names(full_fits), function(model_name) {
        prediction_r2(sampled_predictions$observed, sampled_predictions[[model_name]])
      }, numeric(1))
      heldout_bootstrap[[bootstrap_index]] <- data.table(
        bootstrap = bootstrap_index, model = names(model_r2), heldout_r2 = model_r2,
        heldout_incremental_r2 = model_r2 - model_r2[["covariates_only"]]
      )
    }
    heldout_bootstrap <- rbindlist(heldout_bootstrap)
    heldout_ci <- heldout_bootstrap[, .(
      heldout_r2_q025 = quantile(heldout_r2, 0.025, na.rm = TRUE),
      heldout_r2_q975 = quantile(heldout_r2, 0.975, na.rm = TRUE),
      heldout_incremental_r2_q025 = quantile(heldout_incremental_r2, 0.025, na.rm = TRUE),
      heldout_incremental_r2_q975 = quantile(heldout_incremental_r2, 0.975, na.rm = TRUE)
    ), by = model]
    validation_metrics <- merge(validation_metrics, heldout_ci, by = "model", all.x = TRUE)
    heldout_rows[[paste(trait, validation_name, sep = "__")]] <- validation_metrics

    heldout_wide <- dcast(heldout_bootstrap, bootstrap ~ model, value.var = "heldout_r2")
    heldout_wide[, delta_continuous_minus_four := continuous_three_program - four_group]
    heldout_comparison_rows[[paste(trait, validation_name, sep = "__")]] <- data.table(
      trait = trait, validation_set = validation_name,
      observed_delta_r2_continuous_minus_four =
        validation_metrics[model == "continuous_three_program", heldout_r_squared] -
        validation_metrics[model == "four_group", heldout_r_squared],
      bootstrap_delta_q025 = quantile(heldout_wide$delta_continuous_minus_four, 0.025, na.rm = TRUE),
      bootstrap_delta_q975 = quantile(heldout_wide$delta_continuous_minus_four, 0.975, na.rm = TRUE),
      bootstrap_two_sided_p = 2 * min(
        mean(heldout_wide$delta_continuous_minus_four <= 0, na.rm = TRUE),
        mean(heldout_wide$delta_continuous_minus_four >= 0, na.rm = TRUE)
      )
    )
  }
}

complete_case_flow <- rbindlist(complete_case_rows)
apparent_metrics <- rbindlist(apparent_rows)
bootstrap_summary <- rbindlist(bootstrap_summary_rows)
bootstrap_distribution <- rbindlist(bootstrap_distribution_rows)
cv_distribution <- rbindlist(cv_rows)
cv_fold_distribution <- rbindlist(cv_fold_rows)
cv_comparison <- rbindlist(cv_comparison_rows)
heldout_metrics <- rbindlist(heldout_rows)
heldout_comparison <- rbindlist(heldout_comparison_rows)

cv_summary <- cv_distribution[, .(
  repeated_cv_r2_mean = mean(cv_r_squared, na.rm = TRUE),
  repeated_cv_r2_q025 = quantile(cv_r_squared, 0.025, na.rm = TRUE),
  repeated_cv_r2_q975 = quantile(cv_r_squared, 0.975, na.rm = TRUE),
  repeated_cv_incremental_r2_mean = mean(cv_incremental_r_squared, na.rm = TRUE),
  repeated_cv_incremental_r2_q025 = quantile(cv_incremental_r_squared, 0.025, na.rm = TRUE),
  repeated_cv_incremental_r2_q975 = quantile(cv_incremental_r_squared, 0.975, na.rm = TRUE),
  repeated_cv_rmse_mean = mean(cv_rmse, na.rm = TRUE)
), by = .(trait, model)]

fwrite(complete_case_flow, file.path(output_root, "complete_case_flow.tsv"), sep = "\t")
fwrite(apparent_metrics, file.path(output_root, "apparent_model_metrics.tsv"), sep = "\t")
fwrite(bootstrap_summary, file.path(output_root, "bootstrap_optimism_corrected_metrics.tsv"), sep = "\t")
fwrite(bootstrap_distribution, file.path(output_root, "bootstrap_distribution.tsv.gz"), sep = "\t")
fwrite(cv_summary, file.path(output_root, "repeated_cross_validation_summary.tsv"), sep = "\t")
fwrite(cv_distribution, file.path(output_root, "repeated_cross_validation_distribution.tsv"), sep = "\t")
fwrite(cv_fold_distribution, file.path(output_root, "cross_validation_fold_mse.tsv.gz"), sep = "\t")
fwrite(cv_comparison, file.path(output_root, "formal_cv_continuous_vs_four_group.tsv"), sep = "\t")
fwrite(heldout_metrics, file.path(output_root, "trial_linked_heldout_validation_metrics.tsv"), sep = "\t")
fwrite(heldout_comparison, file.path(output_root, "heldout_continuous_vs_four_group.tsv"), sep = "\t")

decision_summary <- merge(
  bootstrap_summary[model %chin% c("four_group", "continuous_three_program"), .(
    trait, model, adjusted_r_squared, incremental_r_squared,
    incremental_r2_bootstrap_q025, incremental_r2_bootstrap_q975,
    optimism_corrected_incremental_r_squared
  )],
  cv_summary[model %chin% c("four_group", "continuous_three_program"), .(
    trait, model, repeated_cv_incremental_r2_mean,
    repeated_cv_incremental_r2_q025, repeated_cv_incremental_r2_q975
  )],
  by = c("trait", "model")
)
decision_summary <- merge(
  decision_summary,
  heldout_metrics[
    validation_set == "pooled_trial_linked_validation" &
      model %chin% c("four_group", "continuous_three_program"),
    .(trait, model, pooled_heldout_incremental_r2 = heldout_incremental_r_squared,
      pooled_heldout_incremental_q025 = heldout_incremental_r2_q025,
      pooled_heldout_incremental_q975 = heldout_incremental_r2_q975)
  ],
  by = c("trait", "model")
)
fwrite(decision_summary, file.path(output_root, "decision_summary.tsv"), sep = "\t")

capture.output(sessionInfo(), file = file.path(output_root, "sessionInfo.txt"))
message("Clinical landscape validation completed: ", output_root)
