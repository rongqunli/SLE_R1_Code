if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
group_root <- file.path(results_root, "02_group_stratified")
figure_root <- file.path(results_root, "figures")
dir.create(group_root, recursive = TRUE, showWarnings = FALSE)

subtype_map <- unique(group_assignments[, .(subject_id, baseline_subtype = subtype)])
if (anyDuplicated(subtype_map$subject_id)) stop("Baseline grouping map contains duplicated subject IDs")

group_data <- trial_scores[
  analysis_role == "primary" & program %in% primary_programs & time %in% time_levels
]
group_data <- merge(group_data, subtype_map, by = "subject_id", all.x = TRUE)
group_data <- group_data[!is.na(baseline_subtype) & is.finite(score_z)]
group_data[, time := factor(time, levels = time_levels)]
group_data[, treatment := stats::relevel(factor(treatment), ref = cfg$placebo_label)]
group_data[, baseline_subtype := factor(baseline_subtype)]

subtype_audit <- group_data[, .(
  samples = .N,
  subjects = data.table::uniqueN(subject_id),
  baseline_subjects = data.table::uniqueN(subject_id[time == "baseline"]),
  followup_subjects = data.table::uniqueN(subject_id[time != "baseline"])
), by = .(trial, baseline_subtype, treatment, program)]
data.table::fwrite(subtype_audit, file.path(group_root, "subtype_trial_sample_audit.tsv"), sep = "\t")

eligible_group_audit <- group_data[, .(
  subjects = data.table::uniqueN(subject_id),
  placebo_subjects = data.table::uniqueN(subject_id[treatment == cfg$placebo_label]),
  q2w_subjects = data.table::uniqueN(subject_id[treatment == "Tabalumab_Q2W"]),
  q4w_subjects = data.table::uniqueN(subject_id[treatment == "Tabalumab_Q4W"]),
  low_margin_subjects = data.table::uniqueN(subject_id[subject_id %in% group_assignments[low_confidence %in% TRUE, subject_id]])
), by = .(trial, baseline_subtype)]
data.table::fwrite(eligible_group_audit, file.path(group_root, "group_arm_size_audit.tsv"), sep = "\t")

make_did_contrasts <- function(emm, trial, program, placebo, time_levels) {
  grid <- data.table::as.data.table(as.data.frame(emm))
  arms <- setdiff(unique(as.character(grid$treatment)), placebo)
  subtypes <- unique(as.character(grid$baseline_subtype))
  followups <- setdiff(time_levels, "baseline")
  output <- list()
  index <- 0L
  for (subtype_name in subtypes) for (arm in arms) for (followup_time in followups) {
    cells <- list(
      trt_follow = which(grid$treatment == arm & grid$time == followup_time & grid$baseline_subtype == subtype_name),
      trt_base = which(grid$treatment == arm & grid$time == "baseline" & grid$baseline_subtype == subtype_name),
      ctl_follow = which(grid$treatment == placebo & grid$time == followup_time & grid$baseline_subtype == subtype_name),
      ctl_base = which(grid$treatment == placebo & grid$time == "baseline" & grid$baseline_subtype == subtype_name)
    )
    if (!all(lengths(cells) == 1L)) next
    weights <- rep(0, nrow(grid))
    weights[cells$trt_follow] <- 1
    weights[cells$trt_base] <- -1
    weights[cells$ctl_follow] <- -1
    weights[cells$ctl_base] <- 1
    label <- paste(subtype_name, arm, followup_time, sep = "__")
    contrast_table <- as.data.frame(emmeans::contrast(emm, method = setNames(list(weights), label), adjust = "none"))
    statistic_column <- intersect(c("t.ratio", "z.ratio"), names(contrast_table))[1]
    statistic_value <- if (!is.na(statistic_column)) contrast_table[[statistic_column]] else contrast_table$estimate / contrast_table$SE
    df_value <- if ("df" %in% names(contrast_table)) contrast_table$df else Inf
    index <- index + 1L
    output[[index]] <- data.table::data.table(
      trial = trial,
      program = program,
      baseline_subtype = subtype_name,
      treatment = arm,
      time = followup_time,
      estimate = contrast_table$estimate,
      se = contrast_table$SE,
      df = df_value,
      conf_low = contrast_table$estimate - stats::qnorm(0.975) * contrast_table$SE,
      conf_high = contrast_table$estimate + stats::qnorm(0.975) * contrast_table$SE,
      statistic = statistic_value,
      p_value = contrast_table$p.value
    )
  }
  data.table::rbindlist(output, fill = TRUE)
}

fit_group_model <- function(trial_name, program_name, input_data) {
  input_data <- droplevels(input_data)
  subtype_n <- input_data[, data.table::uniqueN(subject_id), by = baseline_subtype]
  keep_subtypes <- subtype_n[V1 >= as.integer(cfg$minimum_subjects_per_group), baseline_subtype]
  input_data <- droplevels(input_data[baseline_subtype %in% keep_subtypes])
  cell_n <- input_data[, .(subjects = data.table::uniqueN(subject_id)), by = .(baseline_subtype, treatment)]
  balanced_subtypes <- cell_n[, .(all_arms_adequate = all(subjects >= as.integer(cfg$minimum_subjects_per_arm))), by = baseline_subtype]
  input_data <- droplevels(input_data[baseline_subtype %in% balanced_subtypes[all_arms_adequate == TRUE, baseline_subtype]])
  complete <- input_data[, .(has_baseline = any(time == "baseline"), has_followup = any(time != "baseline")), by = subject_id]
  complete_ids <- complete[has_baseline & has_followup, subject_id]
  input_data <- droplevels(input_data[subject_id %in% complete_ids])
  if (data.table::uniqueN(input_data$subject_id) < 40L || nlevels(input_data$baseline_subtype) < 2L || nlevels(input_data$treatment) < 2L) return(NULL)

  full_formula <- score_z ~ treatment * time * baseline_subtype + (1 | subject_id)
  reduced_formula <- score_z ~ (treatment + time + baseline_subtype)^2 + (1 | subject_id)
  full <- lme4::lmer(full_formula, data = input_data, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa"))
  reduced <- lme4::lmer(reduced_formula, data = input_data, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa"))
  likelihood_test <- stats::anova(reduced, full)
  model_audit <- data.table::data.table(
    trial = trial_name,
    program = program_name,
    subjects = data.table::uniqueN(input_data$subject_id),
    samples = nrow(input_data),
    subtypes = nlevels(input_data$baseline_subtype),
    treatment_arms = nlevels(input_data$treatment),
    singular = lme4::isSingular(full, tol = 1e-4),
    three_way_chisq = likelihood_test$Chisq[2],
    three_way_p_value = likelihood_test$`Pr(>Chisq)`[2]
  )
  estimated_means <- emmeans::emmeans(full, ~ treatment * time * baseline_subtype, lmer.df = "asymptotic")
  means <- data.table::as.data.table(as.data.frame(estimated_means))
  lower_column <- intersect(c("asymp.LCL", "lower.CL"), names(means))[1]
  upper_column <- intersect(c("asymp.UCL", "upper.CL"), names(means))[1]
  if (is.na(lower_column) || is.na(upper_column)) {
    means[, `:=`(conf_low = emmean - stats::qnorm(0.975) * SE, conf_high = emmean + stats::qnorm(0.975) * SE)]
  } else {
    means[, `:=`(conf_low = get(lower_column), conf_high = get(upper_column))]
  }
  means[, `:=`(trial = trial_name, program = program_name)]
  contrasts <- make_did_contrasts(estimated_means, trial_name, program_name, cfg$placebo_label, time_levels)
  list(model = full, audit = model_audit, means = means, contrasts = contrasts, data = input_data)
}

group_models <- list()
model_index <- 0L
for (trial_name in c("ILLUMINATE-1", "ILLUMINATE-2")) for (program_name in primary_programs) {
  model_index <- model_index + 1L
  group_models[[model_index]] <- fit_group_model(
    trial_name,
    program_name,
    group_data[trial == trial_name & program == program_name]
  )
}
group_models <- Filter(Negate(is.null), group_models)
if (!length(group_models)) stop("No group-stratified model met the sample-size contract")

group_model_audit <- data.table::rbindlist(lapply(group_models, `[[`, "audit"))
group_model_audit[, three_way_fdr := stats::p.adjust(three_way_p_value, method = "BH")]
group_marginal_means <- data.table::rbindlist(lapply(group_models, `[[`, "means"), fill = TRUE)
group_contrasts <- data.table::rbindlist(lapply(group_models, `[[`, "contrasts"), fill = TRUE)
group_contrasts[, fdr := stats::p.adjust(p_value, method = "BH")]
group_contrasts[, direction := data.table::fcase(
  fdr <= as.numeric(cfg$primary_fdr) & estimate < 0, "decrease",
  fdr <= as.numeric(cfg$primary_fdr) & estimate > 0, "increase",
  default = "not_significant"
)]

trial_names <- c("ILLUMINATE-1", "ILLUMINATE-2")
paired_contrasts <- data.table::dcast(
  group_contrasts,
  baseline_subtype + program + treatment + time ~ trial,
  value.var = "estimate"
)
paired_contrasts <- paired_contrasts[is.finite(get(trial_names[1])) & is.finite(get(trial_names[2]))]
cross_trial_consistency <- if (nrow(paired_contrasts)) {
  paired_contrasts[, .(
    shared_contrasts = .N,
    same_direction = sum(sign(get(trial_names[1])) == sign(get(trial_names[2]))),
    same_direction_fraction = mean(sign(get(trial_names[1])) == sign(get(trial_names[2]))),
    spearman_rho = suppressWarnings(stats::cor(get(trial_names[1]), get(trial_names[2]), method = "spearman", use = "complete.obs")),
    mean_effect_trial_1 = mean(get(trial_names[1])),
    mean_effect_trial_2 = mean(get(trial_names[2]))
  ), by = program][, `:=`(trial_1 = trial_names[1], trial_2 = trial_names[2])]
} else data.table::data.table()

data.table::fwrite(group_model_audit, file.path(group_root, "three_way_interaction_tests.tsv"), sep = "\t")
data.table::fwrite(group_marginal_means, file.path(group_root, "model_marginal_means.tsv"), sep = "\t")
data.table::fwrite(group_contrasts, file.path(group_root, "subtype_specific_difference_in_differences.tsv"), sep = "\t")
data.table::fwrite(cross_trial_consistency, file.path(group_root, "cross_trial_direction_consistency.tsv"), sep = "\t")
saveRDS(group_models, file.path(group_root, "trial_specific_mixed_models.rds"))

group_data_audit <- group_data[, .(
  assigned_subjects = data.table::uniqueN(subject_id),
  assigned_samples = .N,
  low_margin_subjects = data.table::uniqueN(subject_id[subject_id %in% group_assignments[low_confidence %in% TRUE, subject_id]])
), by = trial]
data.table::fwrite(group_data_audit, file.path(group_root, "group_model_input_audit.tsv"), sep = "\t")

trajectory_plot <- ggplot2::ggplot(group_marginal_means, ggplot2::aes(time, emmean, color = treatment, group = treatment)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = conf_low, ymax = conf_high, fill = treatment), alpha = 0.10, color = NA) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::facet_grid(trial + baseline_subtype ~ program, scales = "free_y") +
  ggplot2::labs(x = NULL, y = "Estimated program score", color = "Treatment", fill = "Treatment", title = "Grouping-stratified transcriptional pharmacodynamics") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top", strip.text = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(file.path(figure_root, "Figure4C_current_grouping_trajectories.png"), trajectory_plot, width = 15, height = 2.2 + 2.1 * max(1, nlevels(group_marginal_means$baseline_subtype)), dpi = 320, limitsize = FALSE)

effect_plot <- ggplot2::ggplot(group_contrasts, ggplot2::aes(interaction(treatment, time, sep = " | "), interaction(baseline_subtype, program, sep = " | "), fill = estimate)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::geom_point(data = group_contrasts[fdr <= as.numeric(cfg$primary_fdr)], shape = 8, size = 1.6) +
  ggplot2::facet_wrap(~trial, scales = "free_x") +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
  ggplot2::labs(x = "Treatment | follow-up", y = "Baseline grouping | program", fill = "Difference-in-differences", title = "Grouping-specific pharmacodynamic effects") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
ggplot2::ggsave(file.path(figure_root, "Figure4D_current_grouping_effects.png"), effect_plot, width = 14, height = 9, dpi = 320)

cat("Grouping-stratified models and direction checks completed.\n")

