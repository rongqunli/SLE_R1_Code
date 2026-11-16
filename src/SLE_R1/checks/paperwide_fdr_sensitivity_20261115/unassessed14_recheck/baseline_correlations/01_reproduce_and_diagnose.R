if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")


options(digits = 17)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else dirname(normalizePath(sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])))
source_dir <- "F:/program/Hei/SLE_v1/results/reviewer_response/07_bcell_abundance"
samples <- read.delim(file.path(source_dir, "sample_level_b_abundance_and_program_scores.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
old <- read.delim(file.path(source_dir, "b_program_abundance_correlations.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "R_session_info.txt"))
writeLines(deparse(stats:::cor.test.default), file.path(out, "cor_test_default_installed.R"))





as89_diagnostic <- function(n, S, lower = TRUE) {
  q <- round(S) + 2 * lower
  b <- 1 / n
  z <- (6 * (q - 1) * b / (n * n - 1) - 1) * sqrt(n - 1)
  z2 <- z * z
  u <- z * b * (0.2274 + b * (0.2531 + 0.1745 * b) +
      z2 * (-0.0758 + b * (0.1033 + 0.3932 * b) -
      z2 * b * (0.0879 + 0.0151 * b -
      z2 * (0.0072 - 0.0831 * b + z2 * b * (0.0131 - 0.00046 * z2)))))
  correction <- u / exp(z2 / 2)
  normal_tail <- pnorm(z, lower.tail = lower)
  unclipped <- if (lower) normal_tail - correction else normal_tail + correction
  native <- .Call(get("C_pRho", asNamespace("stats")), q, as.integer(n), as.logical(lower))
  data.frame(lower_tail = lower, q_as89 = q, z_as89 = z, normal_tail = normal_tail,
             correction = correction, unclipped_tail = unclipped,
             clipped_tail = max(0, min(1, unclipped)), native_tail = native)
}

results <- list(); diagnostics <- list(); pairs <- list()
for (trial_name in c("ILLUMINATE-1", "ILLUMINATE-2")) {
  baseline <- samples[samples$trial == trial_name & samples$time == "baseline", ]
  stopifnot(!anyNA(baseline$subject_id), !anyDuplicated(baseline$subject_id), !anyDuplicated(baseline$sample_id))
  keep <- is.finite(baseline$official_b_score) & is.finite(baseline$quantiseq_b_fraction_z)
  dat <- baseline[keep, ]
  x <- dat$official_b_score; y <- dat$quantiseq_b_fraction_z; n <- length(x)
  original <- old[old$trial == trial_name & old$subset == "baseline" &
                  old$proxy == "quantiseq_b_fraction_z" & old$method == "spearman", ]
  stopifnot(nrow(original) == 1, n == original$n, all(!duplicated(x)), all(!duplicated(y)))
  default <- cor.test(x, y, method = "spearman")
  explicit_as89 <- cor.test(x, y, method = "spearman", exact = TRUE)
  approximate <- cor.test(x, y, method = "spearman", exact = FALSE)
  rho <- unname(default$estimate)
  t_statistic <- rho / sqrt((1 - rho^2) / (n - 2))
  stable_lower <- 2 * pt(-abs(t_statistic), df = n - 2)
  stable_upper <- 2 * pt(abs(t_statistic), df = n - 2, lower.tail = FALSE)
  stopifnot(abs(rho - original$correlation) < 1e-13,
            default$p.value == original$p_value,
            default$p.value == explicit_as89$p.value,
            abs(log(approximate$p.value) - log(stable_lower)) < 1e-10,
            abs(log(stable_upper) - log(stable_lower)) < 1e-10)
  diagnostic <- rbind(as89_diagnostic(n, unname(default$statistic), TRUE),
                      as89_diagnostic(n, unname(default$statistic), FALSE))
  stopifnot(all(abs(diagnostic$clipped_tail - diagnostic$native_tail) < 1e-14))
  diagnostic$trial <- trial_name
  diagnostics[[trial_name]] <- diagnostic
  pairs[[trial_name]] <- data.frame(trial = trial_name, sample_id = dat$sample_id,
      subject_id = dat$subject_id, official_b_score = x, quantiseq_b_fraction_z = y,
      rank_b_score = rank(x), rank_quantiseq = rank(y))
  results[[trial_name]] <- data.frame(
      hypothesis_id = paste0("baseline_B_fraction_correlation_2::", trial_name,
          "|baseline|quantiseq_b_fraction_z|spearman"),
      trial = trial_name, n_baseline_rows = nrow(baseline), n_complete = n,
      n_excluded_nonfinite = sum(!keep), n_unique_subjects = length(unique(dat$subject_id)),
      n_unique_samples = length(unique(dat$sample_id)),
      tied_x_observations = sum(duplicated(x)), tied_y_observations = sum(duplicated(y)),
      rho_original = original$correlation, rho_reproduced = rho,
      S = unname(default$statistic), p_original = original$p_value,
      p_default_reproduced = default$p.value, p_explicit_as89 = explicit_as89$p.value,
      t_statistic = t_statistic, df = n - 2, p_exact_FALSE = approximate$p.value,
      p_stable_t_lower_tail = stable_lower, p_stable_t_upper_tail = stable_upper,
      p_primary = stable_lower,
      primary_p_method = "two_sided_Spearman_asymptotic_t_exact_FALSE_df_n_minus_2",
      default_zero_cause = if (diagnostic$unclipped_tail[1] < 0) "AS89_Edgeworth_negative_tail_clipped_to_zero" else "other",
      default_normal_tail = diagnostic$normal_tail[1],
      default_correction = diagnostic$correction[1],
      default_unclipped_tail = diagnostic$unclipped_tail[1])
}
write_tsv <- function(x, filename) write.table(x, file.path(out, filename), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write_tsv(do.call(rbind, results), "baseline_correlation_R_recheck.tsv")
write_tsv(do.call(rbind, diagnostics), "AS89_tail_diagnostics.tsv")
write_tsv(do.call(rbind, pairs), "baseline_analysis_pairs.tsv")
print(do.call(rbind, results))
cat("\nAll R reproduction and numerical-tail checks passed.\n")
