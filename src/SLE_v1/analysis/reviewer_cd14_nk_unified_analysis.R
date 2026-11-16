if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

project_root <- "F:/program/Hei/SLE_v1"
output_root <- file.path(project_root, "results", "reviewer_response", "04_cd14_nk_unified")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

score_path <- file.path(project_root, "pipeline", "3.5", "data", "processed", "donor_celltype_program_scores.tsv")
proportion_path <- file.path(project_root, "pipeline", "3.5", "data", "processed", "donor_celltype_proportions.tsv")

scores <- fread(score_path, sep = "\t")
proportions <- fread(proportion_path, sep = "\t")
ifn_program <- "MAGENTA_IFN_ANTIVIRAL"
scores <- scores[
  cohort == "Child" & analysis_set == "primary_ucell" &
    program == ifn_program & condition %in% c("HD", "SLE")
]
proportions <- proportions[
  cohort == "Child" & analysis_set == "primary_ucell" & condition %in% c("HD", "SLE")
]
merge_keys <- c(
  "sample_id", "subject_id", "cohort", "condition", "age", "sex", "batch", "cell_type_analysis"
)
scores <- merge(
  scores,
  proportions[, c(merge_keys, "cell_proportion"), with = FALSE],
  by = merge_keys, all.x = TRUE
)
scores[, condition := factor(condition, levels = c("HD", "SLE"))]
scores[, age_scaled := as.numeric(scale(age))]
scores[, proportion_logit := qlogis(pmin(pmax(cell_proportion, 1e-4), 1 - 1e-4))]
scores[, outcome_z_within_celltype := as.numeric(scale(mean_score)), by = cell_type_analysis]
scores[, outcome_z_global := as.numeric(scale(mean_score))]

coverage <- scores[, .(
  n_HD = uniqueN(subject_id[condition == "HD"]),
  n_SLE = uniqueN(subject_id[condition == "SLE"]),
  median_cells_HD = as.numeric(median(n_cells[condition == "HD"])),
  median_cells_SLE = as.numeric(median(n_cells[condition == "SLE"]))
), by = cell_type_analysis]
eligible <- coverage[n_HD >= 10 & n_SLE >= 20 & median_cells_HD >= 50 & median_cells_SLE >= 50, cell_type_analysis]
all_cell <- copy(scores[cell_type_analysis %in% eligible])
all_cell[, outcome_z := as.numeric(scale(mean_score)), by = cell_type_analysis]
all_cell[, cell_type_analysis := relevel(factor(cell_type_analysis), ref = "NK")]
all_cell <- all_cell[complete.cases(all_cell[, .(
  outcome_z, condition, cell_type_analysis, subject_id, age_scaled, sex, batch, proportion_logit
)])]
all_formula <- outcome_z ~ condition * cell_type_analysis + age_scaled + sex + batch + proportion_logit + (1 | subject_id)
all_fit <- lmer(
  all_formula, data = all_cell, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000), check.conv.singular = "ignore")
)
beta <- fixef(all_fit)
vc <- as.matrix(vcov(all_fit))
terms <- names(beta)
condition_term <- terms[grepl("^conditionSLE$", terms)]
cell_types <- levels(all_cell$cell_type_analysis)
effect_weight <- function(cell_type) {
  w <- setNames(rep(0, length(beta)), terms)
  w[condition_term] <- 1
  if (cell_type != "NK") {
    hit <- terms[grepl(paste0(
      "conditionSLE:cell_type_analysis", cell_type, "$|cell_type_analysis", cell_type, ":conditionSLE$"
    ), terms)]
    if (length(hit) != 1L) stop("Missing all-cell interaction for ", cell_type)
    w[hit] <- 1
  }
  w
}
all_effects <- rbindlist(lapply(cell_types, function(cell_type) {
  w <- effect_weight(cell_type)
  estimate <- sum(w * beta)
  se <- sqrt(as.numeric(t(w) %*% vc %*% w))
  data.table(
    cell_type = cell_type, estimate = estimate, standard_error = se,
    conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se,
    p_value = 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
  )
}))
all_effects[, fdr_cell_effects := p.adjust(p_value, method = "BH")]
fwrite(all_effects, file.path(output_root, "all_celltype_unified_ifn_effects.tsv"), sep = "\t")

fwrite(coverage, file.path(output_root, "all_celltype_coverage.tsv"), sep = "\t")
saveRDS(all_fit, file.path(output_root, "all_celltype_unified_ifn_model.rds"))
