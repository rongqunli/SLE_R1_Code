if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
options(stringsAsFactors = FALSE, width = 140)

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
  library(patchwork)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "analysis", "reviewer_trial_reporting_ancova_permutation.R"), winslash = "/", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "reviewer_response", "08_trial_reporting_ancova_permutation")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, name) {
  fwrite(as.data.table(x), file.path(results_root, name), sep = "\t", na = "NA", quote = FALSE)
}

primary_programs <- c(
  "CYAN_NEUTROPHIL_IMMATURE",
  "GREENYELLOW_B_CELL_COMPLEMENT",
  "MAGENTA_IFN_ANTIVIRAL"
)
coordinate_names <- paste0("coord_", primary_programs)
b_coordinate <- "coord_GREENYELLOW_B_CELL_COMPLEMENT"
coordinate_labels <- c(
  coord_CYAN_NEUTROPHIL_IMMATURE = "Immature-neutrophil",
  coord_GREENYELLOW_B_CELL_COMPLEMENT = "B-cell-enriched",
  coord_MAGENTA_IFN_ANTIVIRAL = "IFN/antiviral"
)
trial_paths <- c(
  `ILLUMINATE-1` = file.path(project_root, "pipeline", "3.4", "data", "processed", "ILLUMINATE1_sample_program_scores.tsv"),
  `ILLUMINATE-2` = file.path(project_root, "pipeline", "3.4", "data", "processed", "ILLUMINATE2_sample_program_scores.tsv")
)
grouping_path <- file.path(project_root, "pipeline", "3.4", "data", "processed", "frozen_grouping_model.rds")
stopifnot(all(file.exists(c(trial_paths, grouping_path))))

grouping_model <- readRDS(grouping_path)
center_values <- grouping_model$center[primary_programs]
scale_values <- grouping_model$scale[primary_programs]
if (any(!is.finite(center_values)) || any(!is.finite(scale_values)) || any(scale_values <= 0)) {
  stop("Frozen coordinate centers/scales are incomplete")
}

standardization_audit <- data.table(
  program = primary_programs,
  coordinate = coordinate_names,
  discovery_center = as.numeric(center_values),
  discovery_scale = as.numeric(scale_values),
  implementation = "(trial score_z - frozen adult-discovery center) / frozen adult-discovery scale"
)
write_tsv(standardization_audit, "coordinate_standardization.tsv")

read_trial <- function(path, trial_name) {
  score_long <- fread(path, sep = "\t")
  score_long <- score_long[
    analysis_role == "primary" & program %in% primary_programs &
      time %in% c("baseline", "week16", "week52") & is.finite(score_z)
  ]
  wide <- dcast(
    score_long,
    sample_id + subject_id + time + treatment ~ program,
    value.var = "score_z"
  )
  if (anyDuplicated(wide[, .(subject_id, time)])) stop("Duplicate subject-time rows in ", trial_name)
  raw <- as.matrix(wide[, ..primary_programs])
  coords <- sweep(sweep(raw, 2, center_values, `-`), 2, scale_values, `/`)
  colnames(coords) <- coordinate_names
  wide[, (coordinate_names) := as.data.table(coords)]
  wide[, trial := trial_name]
  wide
}

trial_wide <- rbindlist(Map(read_trial, trial_paths, names(trial_paths)), use.names = TRUE, fill = TRUE)
setcolorder(trial_wide, c("trial", "sample_id", "subject_id", "time", "treatment", primary_programs, coordinate_names))

arm_consistency <- trial_wide[, .(n_treatments = uniqueN(treatment)), by = .(trial, subject_id)]
if (any(arm_consistency$n_treatments != 1L)) stop("Treatment arm is inconsistent within participant")

arm_subjects <- unique(trial_wide[, .(trial, subject_id, treatment)])
availability <- unique(trial_wide[, .(trial, subject_id, treatment, time)])
flow_rows <- list()
for (tr in names(trial_paths)) {
  for (arm in c("Placebo", "Tabalumab_Q2W", "Tabalumab_Q4W")) {
    subjects <- arm_subjects[trial == tr & treatment == arm, subject_id]
    base_ids <- availability[trial == tr & treatment == arm & time == "baseline", subject_id]
    for (tt in c("week16", "week52")) {
      follow_ids <- availability[trial == tr & treatment == arm & time == tt, subject_id]
      flow_rows[[length(flow_rows) + 1L]] <- data.table(
        trial = tr,
        treatment = arm,
        followup = tt,
        transcriptomic_cohort_subjects = length(subjects),
        baseline_available = length(base_ids),
        followup_available = length(follow_ids),
        paired_baseline_followup = length(intersect(base_ids, follow_ids)),
        baseline_missing = length(setdiff(subjects, base_ids)),
        followup_missing_among_baseline = length(setdiff(base_ids, follow_ids))
      )
    }
  }
}
trial_flow <- rbindlist(flow_rows)
write_tsv(trial_flow, "trial_transcriptomic_flow.tsv")




parent_trial_randomization <- data.table(
  trial = rep(c("ILLUMINATE-1", "ILLUMINATE-2"), each = 3L),
  treatment = rep(c("Placebo", "Tabalumab_Q2W", "Tabalumab_Q4W"), 2L),
  parent_trial_randomized = c(388L, 387L, 389L, 376L, 372L, 376L),
  trial_total_randomized = rep(c(1164L, 1124L), each = 3L),
  trial_registration = rep(c("NCT01196091", "NCT01205438"), each = 3L),
  primary_report_doi = rep(c("10.1136/annrheumdis-2015-207653", "10.1136/annrheumdis-2015-207654"), each = 3L),
  background_therapy = "study drug or placebo added to continued standard-of-care SLE therapy"
)
write_tsv(parent_trial_randomization, "parent_trial_randomization.tsv")
trial_flow_with_parent <- merge(trial_flow, parent_trial_randomization, by = c("trial", "treatment"), all.x = TRUE)
trial_flow_with_parent[, molecular_subset_fraction_of_randomized := transcriptomic_cohort_subjects / parent_trial_randomized]
write_tsv(trial_flow_with_parent, "trial_flow_parent_and_molecular.tsv")

all_available_audit <- trial_wide[, .(
  n_samples = .N,
  n_subjects = uniqueN(subject_id)
), by = .(trial, treatment, time)][order(trial, treatment, match(time, c("baseline", "week16", "week52")))]
write_tsv(all_available_audit, "all_available_mixed_model_audit.tsv")

metadata_availability <- data.table(
  requested_item = c(
    "molecular sample and visit availability", "assigned arm", "withdrawal/dropout reason",
    "concomitant medication", "clinical baseline covariates", "parent-trial randomized population"
  ),
  available_in_public_GEO_molecular_metadata = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  handling = c(
    "reported by trial, arm, and visit", "reported from propagated subject-level assignment",
    "not supplied; no reason-specific attribution", "not supplied; cannot adjust in molecular subset",
    "not supplied; baseline balance limited to frozen molecular coordinates",
    "reported from the primary trial publications and kept separate from GEO molecular-analysis counts"
  )
)
write_tsv(metadata_availability, "metadata_availability_and_handling.tsv")

analysis_set_definitions <- data.table(
  analysis = c(
    "longitudinal mixed model", "change-score vector sensitivity",
    "follow-up ANCOVA", "baseline-by-treatment ANCOVA", "treatment-label permutation"
  ),
  analysis_set = c(
    "all available molecular samples across baseline, week 16, and week 52",
    "complete baseline/follow-up pairs at the specified visit",
    "complete baseline/follow-up pairs at the specified visit",
    "complete baseline/follow-up pairs at the specified visit",
    "same complete-pair contrast; arm labels permuted with arm sizes fixed"
  ),
  missing_data_handling = c(
    "likelihood-based estimation; no expression-value imputation",
    "no imputation", "no imputation", "no imputation", "not applicable"
  ),
  estimand = c(
    "treatment-by-visit difference from placebo",
    "drug-minus-placebo difference in mean within-person change",
    "placebo-adjusted follow-up mean conditional on baseline coordinate",
    "placebo-adjusted treatment-by-baseline-coordinate interaction",
    "randomization-null distribution of vector selectivity"
  )
)
write_tsv(analysis_set_definitions, "analysis_set_definitions.tsv")

baseline_rows <- trial_wide[time == "baseline"]
baseline_summary <- melt(
  baseline_rows,
  id.vars = c("trial", "subject_id", "treatment"),
  measure.vars = coordinate_names,
  variable.name = "coordinate",
  value.name = "value"
)[, .(n = .N, mean = mean(value), sd = sd(value), median = median(value), q1 = quantile(value, 0.25), q3 = quantile(value, 0.75)),
  by = .(trial, treatment, coordinate)]
write_tsv(baseline_summary, "baseline_coordinate_balance_by_arm.tsv")

pairwise_balance <- list()
for (tr in names(trial_paths)) {
  for (coord in coordinate_names) {
    dat <- baseline_rows[trial == tr, .(treatment, value = get(coord))]
    for (dose in c("Tabalumab_Q2W", "Tabalumab_Q4W")) {
      a <- dat[treatment == dose, value]
      p <- dat[treatment == "Placebo", value]
      pooled_sd <- sqrt(((length(a) - 1) * var(a) + (length(p) - 1) * var(p)) / (length(a) + length(p) - 2))
      pairwise_balance[[length(pairwise_balance) + 1L]] <- data.table(
        trial = tr, dose = dose, coordinate = coord,
        n_dose = length(a), n_placebo = length(p),
        standardized_mean_difference = (mean(a) - mean(p)) / pooled_sd
      )
    }
  }
}
pairwise_balance <- rbindlist(pairwise_balance)
write_tsv(pairwise_balance, "baseline_coordinate_pairwise_smd.tsv")

make_paired <- function(trial_name, followup_time) {
  base <- copy(trial_wide[trial == trial_name & time == "baseline", c("subject_id", "treatment", coordinate_names), with = FALSE])
  follow <- copy(trial_wide[trial == trial_name & time == followup_time, c("subject_id", "treatment", coordinate_names), with = FALSE])
  setnames(base, coordinate_names, paste0(coordinate_names, "_baseline"))
  setnames(follow, coordinate_names, paste0(coordinate_names, "_followup"))
  joined <- merge(base, follow, by = c("subject_id", "treatment"), all = FALSE)
  for (coord in coordinate_names) {
    joined[, (paste0(coord, "_change")) := get(paste0(coord, "_followup")) - get(paste0(coord, "_baseline"))]
  }
  joined[, `:=`(trial = trial_name, time = followup_time)]
  joined
}

paired_data <- rbindlist(lapply(names(trial_paths), function(tr) {
  rbindlist(lapply(c("week16", "week52"), function(tt) make_paired(tr, tt)), use.names = TRUE)
}), use.names = TRUE)
write_tsv(paired_data, "complete_pair_coordinate_data.tsv")

extract_hc3 <- function(fit, term) {
  ct <- coeftest(fit, vcov. = vcovHC(fit, type = "HC3"))
  if (!term %in% rownames(ct)) return(data.table(estimate = NA_real_, std_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_))
  est <- unname(ct[term, "Estimate"])
  se <- unname(ct[term, "Std. Error"])
  data.table(
    estimate = est,
    std_error = se,
    conf_low = est - qt(0.975, df.residual(fit)) * se,
    conf_high = est + qt(0.975, df.residual(fit)) * se,
    p_value = unname(ct[term, "Pr(>|t|)"])
  )
}

ancova_rows <- list()
moderation_rows <- list()
for (tr in names(trial_paths)) {
  for (tt in c("week16", "week52")) {
    for (dose in c("Tabalumab_Q2W", "Tabalumab_Q4W")) {
      contrast <- copy(paired_data[trial == tr & time == tt & treatment %in% c("Placebo", dose)])
      contrast[, treatment_indicator := as.integer(treatment == dose)]
      n_dose <- sum(contrast$treatment_indicator == 1L)
      n_placebo <- sum(contrast$treatment_indicator == 0L)
      for (coord in coordinate_names) {
        dat <- contrast[, .(
          followup_value = get(paste0(coord, "_followup")),
          baseline_value = get(paste0(coord, "_baseline")),
          treatment_indicator
        )]
        fit <- lm(followup_value ~ baseline_value + treatment_indicator, data = dat)
        row <- extract_hc3(fit, "treatment_indicator")
        row[, `:=`(
          trial = tr, time = tt, dose = dose, coordinate = coord,
          n_dose = n_dose, n_placebo = n_placebo,
          baseline_coefficient = unname(coef(fit)["baseline_value"]),
          method = "follow-up ANCOVA with HC3 covariance"
        )]
        ancova_rows[[length(ancova_rows) + 1L]] <- row
      }

      dat_b <- contrast[, .(
        followup_value = get(paste0(b_coordinate, "_followup")),
        baseline_value = get(paste0(b_coordinate, "_baseline")),
        treatment_indicator
      )]
      dat_b[, baseline_value_z := as.numeric(scale(baseline_value))]
      fit_int <- lm(followup_value ~ baseline_value_z * treatment_indicator, data = dat_b)
      int_row <- extract_hc3(fit_int, "baseline_value_z:treatment_indicator")
      int_row[, `:=`(
        trial = tr, time = tt, dose = dose,
        n_dose = n_dose, n_placebo = n_placebo,
        baseline_reference = "within-contrast mean; coefficient per 1 SD higher baseline B-cell coordinate",
        interpretation = "additional tabalumab-versus-placebo follow-up difference per baseline SD"
      )]
      moderation_rows[[length(moderation_rows) + 1L]] <- int_row
    }
  }
}
ancova_results <- rbindlist(ancova_rows, fill = TRUE)
ancova_results[, fdr := p.adjust(p_value, method = "BH")]
setcolorder(ancova_results, c("trial", "time", "dose", "coordinate", "estimate", "std_error", "conf_low", "conf_high", "p_value", "fdr"))
write_tsv(ancova_results, "ancova_coordinate_effects_trial_specific.tsv")

moderation_results <- rbindlist(moderation_rows, fill = TRUE)
moderation_results[, fdr := p.adjust(p_value, method = "BH")]
setcolorder(moderation_results, c("trial", "time", "dose", "estimate", "std_error", "conf_low", "conf_high", "p_value", "fdr"))
write_tsv(moderation_results, "ancova_bcell_baseline_interactions_trial_specific.tsv")

fixed_effect_meta <- function(dat, group_cols) {
  dat[, {
    ok <- is.finite(estimate) & is.finite(std_error) & std_error > 0
    w <- 1 / std_error[ok]^2
    est <- sum(w * estimate[ok]) / sum(w)
    se <- sqrt(1 / sum(w))
    q <- sum(w * (estimate[ok] - est)^2)
    data.table(
      estimate = est, std_error = se,
      conf_low = est - qnorm(0.975) * se,
      conf_high = est + qnorm(0.975) * se,
      p_value = 2 * pnorm(-abs(est / se)),
      heterogeneity_q = q,
      heterogeneity_df = sum(ok) - 1L,
      heterogeneity_p = if (sum(ok) > 1L) pchisq(q, df = sum(ok) - 1L, lower.tail = FALSE) else NA_real_,
      n_trials = sum(ok)
    )
  }, by = group_cols]
}

ancova_meta <- fixed_effect_meta(ancova_results, c("time", "dose", "coordinate"))
ancova_meta[, fdr := p.adjust(p_value, method = "BH")]
write_tsv(ancova_meta, "ancova_coordinate_effects_fixed_effect_meta.tsv")

moderation_meta <- fixed_effect_meta(moderation_results, c("time", "dose"))
moderation_meta[, fdr := p.adjust(p_value, method = "BH")]
write_tsv(moderation_meta, "ancova_bcell_baseline_interactions_fixed_effect_meta.tsv")

ancova_vector <- function(dat) {
  out <- setNames(numeric(length(coordinate_names)), coordinate_names)
  for (coord in coordinate_names) {
    y <- dat[[paste0(coord, "_followup")]]
    x <- dat[[paste0(coord, "_baseline")]]
    z <- dat$treatment_indicator
    out[[coord]] <- unname(lm.fit(cbind(1, x, z), y)$coefficients[[3]])
  }
  out
}

change_vector <- function(dat) {
  out <- vapply(coordinate_names, function(coord) {
    values <- dat[[paste0(coord, "_change")]]
    mean(values[dat$treatment_indicator == 1L]) - mean(values[dat$treatment_indicator == 0L])
  }, numeric(1))
  setNames(out, coordinate_names)
}

n_boot <- 3000L
n_perm <- 5000L
set.seed(20260811L)
component_rows <- list()
ancova_bootstrap_draws <- list()
for (tr in names(trial_paths)) {
  for (tt in c("week16", "week52")) {
    for (dose in c("Tabalumab_Q2W", "Tabalumab_Q4W")) {
      dat <- copy(paired_data[trial == tr & time == tt & treatment %in% c("Placebo", dose)])
      dat[, treatment_indicator := as.integer(treatment == dose)]
      drug_idx <- which(dat$treatment_indicator == 1L)
      placebo_idx <- which(dat$treatment_indicator == 0L)
      n_dose <- length(drug_idx)
      n_placebo <- length(placebo_idx)
      for (method in c("ANCOVA", "change_score")) {
        vector_fun <- if (method == "ANCOVA") ancova_vector else change_vector
        point <- vector_fun(dat)
        boot_vectors <- matrix(NA_real_, nrow = n_boot, ncol = length(coordinate_names), dimnames = list(NULL, coordinate_names))
        for (b in seq_len(n_boot)) {
          sampled <- rbindlist(list(dat[sample(placebo_idx, n_placebo, replace = TRUE)],dat[sample(drug_idx, n_dose, replace = TRUE)]))
          boot_vectors[b, ] <- vector_fun(sampled)
        }
        if (method == "ANCOVA") ancova_bootstrap_draws[[paste(tr, tt, dose, sep = "__")]] <- boot_vectors
        for (p in seq_len(n_perm)) invisible(sample(dat$treatment_indicator, replace = FALSE))
        for (coord in coordinate_names) {
          component_rows[[length(component_rows) + 1L]] <- data.table(trial = tr, time = tt, dose = dose, method = method, coordinate = coord, estimate = point[[coord]], conf_low = quantile(boot_vectors[, coord], 0.025, na.rm = TRUE), conf_high = quantile(boot_vectors[, coord], 0.975, na.rm = TRUE), probability_negative = mean(boot_vectors[, coord] < 0, na.rm = TRUE), n_dose = n_dose, n_placebo = n_placebo)
        }
      }
    }
  }
}
write_tsv(rbindlist(component_rows), "vector_components_bootstrap.tsv")
saveRDS(ancova_bootstrap_draws, file.path(results_root, "ancova_vector_bootstrap_draws.rds"), compress = TRUE)
