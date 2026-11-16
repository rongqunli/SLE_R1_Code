if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
overall_root <- file.path(results_root, "01_overall")
figure_root <- file.path(results_root, "figures")
dir.create(overall_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

extract_coefficient <- function(model, term, program_name, model_family, n_subjects, n_samples) {
  coefficients <- coef(summary(model))
  result <- data.table::data.table(
    program = program_name,
    model_family = model_family,
    term = term,
    estimate = NA_real_,
    std_error = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    n_subjects = n_subjects,
    n_samples = n_samples,
    p_method = "two-sided normal approximation"
  )
  if (!(term %in% rownames(coefficients))) return(result)
  estimate_value <- unname(coefficients[term, "Estimate"])
  std_error_value <- unname(coefficients[term, "Std. Error"])
  statistic_value <- estimate_value / std_error_value
  result[, `:=`(
    estimate = estimate_value,
    std_error = std_error_value,
    conf_low = estimate_value - 1.96 * std_error_value,
    conf_high = estimate_value + 1.96 * std_error_value,
    statistic = statistic_value,
    p_value = 2 * stats::pnorm(abs(statistic_value), lower.tail = FALSE)
  )]
  result
}

overall_data <- trial_scores[
  analysis_role == "primary" & program %in% primary_programs & time %in% time_levels & is.finite(score_z)
]
overall_data[, time := stats::relevel(factor(time, levels = time_levels), ref = "baseline")]
overall_data[, treatment := stats::relevel(factor(treatment), ref = cfg$placebo_label)]

target_terms <- c(
  "timeweek16",
  "timeweek52",
  "timeweek16:treatmentTabalumab_Q2W",
  "timeweek52:treatmentTabalumab_Q2W",
  "timeweek16:treatmentTabalumab_Q4W",
  "timeweek52:treatmentTabalumab_Q4W"
)

overall_models <- list()
overall_results <- list()
paired_results <- list()
result_index <- 0L
paired_index <- 0L

for (trial_name in unique(overall_data$trial)) {
  trial_dir <- file.path(overall_root, trial_name)
  dir.create(trial_dir, recursive = TRUE, showWarnings = FALSE)
  trial_models <- list()
  trial_result <- data.table::data.table()
  trial_paired <- data.table::data.table()

  for (program_name in unique(overall_data$program)) {
    dat <- droplevels(overall_data[trial == trial_name & program == program_name])
    fit <- lme4::lmer(
      score_z ~ time * treatment + (1 | subject_id),
      data = dat,
      REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 200000),
        check.conv.singular = "ignore"
      )
    )
    trial_models[[program_name]] <- fit
    rows <- data.table::rbindlist(lapply(target_terms, function(term) {
      extract_coefficient(
        fit, term, program_name, "longitudinal_mixed",
        data.table::uniqueN(dat$subject_id), nrow(dat)
      )
    }))
    rows[, analysis_role := "primary"]
    trial_result <- data.table::rbindlist(list(trial_result, rows), use.names = TRUE)

    baseline <- dat[time == "baseline", .(subject_id, baseline_score = score_z)]
    for (target_time in c("week16", "week52")) {
      followup <- dat[time == target_time, .(subject_id, treatment, followup_score = score_z)]
      paired <- merge(followup, baseline, by = "subject_id")
      paired[, delta := followup_score - baseline_score]
      paired[, treatment := stats::relevel(factor(treatment), ref = cfg$placebo_label)]
      paired_fit <- stats::lm(delta ~ treatment, data = paired)
      paired_rows <- data.table::rbindlist(lapply(
        c("treatmentTabalumab_Q2W", "treatmentTabalumab_Q4W"),
        function(term) extract_coefficient(
          paired_fit, term, program_name, paste0("paired_change_", target_time),
          data.table::uniqueN(paired$subject_id), nrow(paired)
        )
      ))
      paired_rows[, analysis_role := "primary"]
      trial_paired <- data.table::rbindlist(list(trial_paired, paired_rows), use.names = TRUE)
    }
  }

  trial_result[, fdr := stats::p.adjust(p_value, method = "BH"), by = .(term, analysis_role)]
  trial_paired[, fdr := stats::p.adjust(p_value, method = "BH"), by = .(model_family, analysis_role)]
  data.table::fwrite(trial_result, file.path(trial_dir, "mixed_model_results.tsv"), sep = "\t")
  data.table::fwrite(trial_paired, file.path(trial_dir, "paired_change_sensitivity.tsv"), sep = "\t")

  overall_models[[trial_name]] <- trial_models
  result_index <- result_index + 1L
  overall_results[[result_index]] <- data.table::copy(trial_result)[, trial := trial_name]
  paired_index <- paired_index + 1L
  paired_results[[paired_index]] <- data.table::copy(trial_paired)[, trial := trial_name]
}

overall_coefficients <- data.table::rbindlist(overall_results, use.names = TRUE, fill = TRUE)
paired_change_results <- data.table::rbindlist(paired_results, use.names = TRUE, fill = TRUE)
data.table::fwrite(overall_coefficients, file.path(overall_root, "trial_mixed_model_coefficients.tsv"), sep = "\t")
data.table::fwrite(paired_change_results, file.path(overall_root, "paired_change_summary.tsv"), sep = "\t")
saveRDS(overall_models, file.path(overall_root, "overall_trial_mixed_models.rds"))

bcell_program <- "GREENYELLOW_B_CELL_COMPLEMENT"
trajectory_summary <- overall_data[program == bcell_program, .(
  n_subjects = data.table::uniqueN(subject_id),
  n_samples = .N,
  mean_score = mean(score_z),
  standard_error = stats::sd(score_z) / sqrt(.N),
  conf_low = mean(score_z) - stats::qnorm(0.975) * stats::sd(score_z) / sqrt(.N),
  conf_high = mean(score_z) + stats::qnorm(0.975) * stats::sd(score_z) / sqrt(.N)
), by = .(trial, treatment, time)]
trajectory_summary[, time := factor(time, levels = time_levels)]
data.table::setorder(trajectory_summary, trial, treatment, time)
data.table::fwrite(trajectory_summary, file.path(overall_root, "overall_trajectory_summary.tsv"), sep = "\t")

treatment_colors <- c(Placebo = "#6C757D", Tabalumab_Q2W = "#1B7F79", Tabalumab_Q4W = "#C65D3B")
trajectory_plot <- ggplot2::ggplot(
  trajectory_summary,
  ggplot2::aes(x = time, y = mean_score, group = treatment, color = treatment)
) +
  ggplot2::geom_hline(yintercept = 0, color = "#BFC3C7", linewidth = 0.4) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.08, linewidth = 0.55) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::facet_wrap(~trial, nrow = 1) +
  ggplot2::scale_color_manual(values = treatment_colors) +
  ggplot2::labs(
    title = "B-cell/complement program trajectories in two randomized trials",
    subtitle = "Observed means and 95% confidence intervals; participants are the longitudinal analysis unit",
    x = NULL,
    y = "Mean standardized program score",
    color = "Treatment"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(
  file.path(figure_root, "Figure4A_current_overall_bcell_trajectories.png"),
  trajectory_plot,
  width = 9,
  height = 4.8,
  dpi = 320
)

cat("Overall trial mixed models and paired-change sensitivities completed.\n")
