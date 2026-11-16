if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
cross_root <- file.path(results_root, "03_cross_trial")
figure_root <- file.path(results_root, "figures")
dir.create(cross_root, recursive = TRUE, showWarnings = FALSE)

trial_models <- data.table::rbindlist(lapply(
  c("ILLUMINATE-1", "ILLUMINATE-2"),
  function(trial_name) {
    x <- data.table::fread(file.path(results_root, "01_overall", trial_name, "mixed_model_results.tsv"), sep = "\t")
    x[, trial := trial_name]
    x
  }
), use.names = TRUE, fill = TRUE)
trial_models <- trial_models[grepl("^timeweek(16|52):treatmentTabalumab_(Q2W|Q4W)$", term)]
trial_models[, week := sub("^timeweek([0-9]+):.*$", "\\1", term)]
trial_models[, dose := sub("^.*treatmentTabalumab_", "", term)]

meta_summary <- trial_models[, {
  ordering <- order(match(trial, c("ILLUMINATE-1", "ILLUMINATE-2")))
  estimate_ordered <- estimate[ordering]
  se_ordered <- std_error[ordering]
  trial_ordered <- trial[ordering]
  weight <- 1 / se_ordered^2
  pooled <- sum(weight * estimate_ordered) / sum(weight)
  pooled_se <- sqrt(1 / sum(weight))
  q_statistic <- sum(weight * (estimate_ordered - pooled)^2)
  q_df <- length(estimate_ordered) - 1L
  data.table::data.table(
    n_trials = length(estimate_ordered),
    trial_1 = trial_ordered[1],
    trial_1_estimate = estimate_ordered[1],
    trial_1_standard_error = se_ordered[1],
    trial_2 = trial_ordered[2],
    trial_2_estimate = estimate_ordered[2],
    trial_2_standard_error = se_ordered[2],
    meta_estimate = pooled,
    meta_standard_error = pooled_se,
    conf_low = pooled - stats::qnorm(0.975) * pooled_se,
    conf_high = pooled + stats::qnorm(0.975) * pooled_se,
    z_statistic = pooled / pooled_se,
    p_value = 2 * stats::pnorm(abs(pooled / pooled_se), lower.tail = FALSE),
    q_statistic = q_statistic,
    q_df = q_df,
    q_p_value = stats::pchisq(q_statistic, df = q_df, lower.tail = FALSE),
    i2_percent = if (q_statistic > 0) max(0, (q_statistic - q_df) / q_statistic) * 100 else 0
  )
}, by = .(program, dose, week)]
meta_summary[, dose_order := match(dose, c("Q2W", "Q4W"))]
meta_summary[, week_order := as.integer(week)]
data.table::setorder(meta_summary, program, dose_order, week_order)
meta_summary[, c("dose_order", "week_order") := NULL]
meta_summary[, fdr := stats::p.adjust(p_value, method = "BH")]
meta_summary <- meta_summary[, .(
  program, dose, week, n_trials, trial_1, trial_1_estimate, trial_1_standard_error,
  trial_2, trial_2_estimate, trial_2_standard_error, meta_estimate, meta_standard_error,
  conf_low, conf_high, z_statistic, p_value, q_statistic, q_df, q_p_value, i2_percent, fdr
)]
data.table::fwrite(meta_summary, file.path(cross_root, "cross_trial_fixed_effect_meta.tsv"), sep = "\t")

legacy_meta <- trial_models[, {
  ordering <- order(match(trial, c("ILLUMINATE-1", "ILLUMINATE-2")))
  estimate_ordered <- estimate[ordering]
  se_ordered <- std_error[ordering]
  trial_ordered <- c("GSE88885", "GSE88886")
  weight <- 1 / se_ordered^2
  pooled <- sum(weight * estimate_ordered) / sum(weight)
  pooled_se <- sqrt(1 / sum(weight))
  q_statistic <- sum(weight * (estimate_ordered - pooled)^2)
  data.table::data.table(
    k = length(estimate_ordered),
    pooled_estimate = pooled,
    pooled_std_error = pooled_se,
    conf_low = pooled - 1.96 * pooled_se,
    conf_high = pooled + 1.96 * pooled_se,
    z_value = pooled / pooled_se,
    p_value = 2 * stats::pnorm(abs(pooled / pooled_se), lower.tail = FALSE),
    heterogeneity_q = q_statistic,
    heterogeneity_df = length(estimate_ordered) - 1L,
    i2_percent = if (q_statistic > 0) max(0, (q_statistic - 1) / q_statistic) * 100 else 0,
    trial_estimates = paste(paste0(trial_ordered, "=", sprintf("%.4f", estimate_ordered)), collapse = ";")
  )
}, by = .(program, term)]
legacy_meta[, analysis_role := "primary"]
legacy_meta[, fdr := stats::p.adjust(p_value, method = "BH")]
legacy_meta <- legacy_meta[, .(
  program, analysis_role, term, k, pooled_estimate, pooled_std_error, conf_low, conf_high,
  z_value, p_value, heterogeneity_q, heterogeneity_df, i2_percent, trial_estimates, fdr
)]
data.table::setorder(legacy_meta, program, term)
data.table::fwrite(legacy_meta, file.path(cross_root, "cross_trial_pharmacodynamic_meta.tsv"), sep = "\t")

bcell_meta <- meta_summary[program == "GREENYELLOW_B_CELL_COMPLEMENT"]
bcell_meta[, contrast := paste0(dose, ", week ", week)]
bcell_meta[, contrast := factor(contrast, levels = rev(c("Q2W, week 16", "Q4W, week 16", "Q2W, week 52", "Q4W, week 52")))]
forest_plot <- ggplot2::ggplot(bcell_meta, ggplot2::aes(x = meta_estimate, y = contrast)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#7B7D7D") +
  ggplot2::geom_segment(ggplot2::aes(x = conf_low, xend = conf_high, yend = contrast), color = "#2F7D67", linewidth = 1.2) +
  ggplot2::geom_point(color = "#2F7D67", size = 3) +
  ggplot2::geom_text(ggplot2::aes(label = paste0("FDR=", format.pval(fdr, digits = 3, eps = 1e-3), "; I2=", sprintf("%.0f%%", i2_percent))), x = -0.03, hjust = 1, size = 3.1, color = "#303030") +
  ggplot2::scale_x_continuous(limits = c(-0.72, 0.02), breaks = seq(-0.7, 0, 0.1)) +
  ggplot2::labs(title = "Cross-trial fixed-effect synthesis", subtitle = "Treatment-by-time effects on the B-cell/complement program; bars are 95% confidence intervals", x = "Standardized treatment-by-time effect", y = NULL) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(figure_root, "Figure4B_current_cross_trial_meta_forest.png"), forest_plot, width = 8, height = 4.8, dpi = 320)

comparison_families <- data.table::data.table(
  family = c("overall_meta", "global_three_way", "eligible_group_contrasts", "shared_group_direction_checks"),
  definition = c(
    "Three programs x two doses x two follow-up times; 12 prespecified overall comparisons; BH FDR across all 12.",
    "Three programs x two trials; six treatment-by-time-by-group likelihood-ratio tests; BH FDR across all six.",
    "All estimable group-specific dose-by-time difference-in-differences in both trials; BH FDR across the eligible contrast family.",
    "Eight B-cell/complement contrasts shared by IMMUNE_LOW and IFN_HIGH across two trials; direction agreement is descriptive and is not a second meta-analysis."
  ),
  interpretation = c("Overall molecular pharmacodynamics", "Group-heterogeneity screen", "Supportive grouping-dependent amplitude", "Directional concordance only")
)
data.table::fwrite(comparison_families, file.path(cross_root, "comparison_families.tsv"), sep = "\t")
cat("Cross-trial fixed-effect synthesis completed.\n")
