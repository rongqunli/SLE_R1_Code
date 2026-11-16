if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- if (length(file_arg)) normalizePath(file_arg[[1]], winslash = "/", mustWork = TRUE) else normalizePath(file.path(getwd(), "analysis", "followup_landscape_mechanisms.R"), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "exploratory_program_landscape")
followup_root <- file.path(results_root, "04_followup_mechanisms")
for (directory in c("01_bcell_localization", "02_treatment_vectors", "03_ifn_cell_mechanism", "figures", "00_audit")) dir.create(file.path(followup_root, directory), recursive = TRUE, showWarnings = FALSE)
required_packages <- c("data.table", "ggplot2", "lme4")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
options(stringsAsFactors = FALSE, width = 180)
primary_programs <- c("CYAN_NEUTROPHIL_IMMATURE", "GREENYELLOW_B_CELL_COMPLEMENT", "MAGENTA_IFN_ANTIVIRAL")
bcell_program <- "GREENYELLOW_B_CELL_COMPLEMENT"
ifn_program <- "MAGENTA_IFN_ANTIVIRAL"
low_group <- "S02_IMMUNE_LOW"
input_paths <- c(
  donor_scores = file.path(project_root, "pipeline", "3.5", "data", "processed", "donor_celltype_program_scores.tsv"),
  donor_proportions = file.path(project_root, "pipeline", "3.5", "data", "processed", "donor_celltype_proportions.tsv"),
  program_matching = file.path(project_root, "pipeline", "3.5", "data", "processed", "frozen_program_gene_matching.tsv"),
  pseudobulk = file.path(project_root, "pipeline", "3.5", "data", "processed", "pseudobulk_all_genes.rds")
)
if (any(!file.exists(input_paths))) stop("Missing inputs: ", paste(names(input_paths)[!file.exists(input_paths)], collapse = ", "))
safe_numeric <- function(x) suppressWarnings(as.numeric(x))
available_covariates <- function(dat, candidates) candidates[vapply(candidates, function(v) length(unique(dat[[v]][!is.na(dat[[v]])])) > 1L, logical(1))]
formula_from_terms <- function(outcome, terms) if (length(terms)) stats::as.formula(paste(outcome, "~", paste(terms, collapse = " + "))) else stats::as.formula(paste(outcome, "~ 1"))
write_tsv <- function(x, relative_path) data.table::fwrite(x, file.path(followup_root, relative_path), sep = "\t", na = "NA")

donor_scores <- data.table::fread(input_paths[["donor_scores"]], sep = "\t")
donor_proportions <- data.table::fread(input_paths[["donor_proportions"]], sep = "\t")
donor_scores <- donor_scores[analysis_set == "primary_ucell"]
donor_scores[, mean_score := safe_numeric(mean_score)]
donor_scores[, condition := factor(condition, levels = c("HD", "SLE"))]
donor_scores[, cohort := factor(cohort, levels = c("Child", "Adult"))]
donor_scores[, batch := factor(batch)]
donor_scores[, sex := factor(sex)]
donor_proportions <- donor_proportions[analysis_set == "primary_ucell", .(sample_id, cell_type_analysis, cell_proportion = safe_numeric(cell_proportion))]
donor_scores <- merge(donor_scores, donor_proportions, by = c("sample_id", "cell_type_analysis"), all.x = TRUE)

bcell_scores <- donor_scores[program == bcell_program]
bcell_cell_summary <- bcell_scores[, .(n_profiles = .N, n_donors = data.table::uniqueN(subject_id), mean_score = mean(mean_score, na.rm = TRUE), standard_deviation = stats::sd(mean_score, na.rm = TRUE), mean_cell_proportion = mean(cell_proportion, na.rm = TRUE)), by = .(cohort, condition, cell_type_analysis)]
write_tsv(bcell_cell_summary, "01_bcell_localization/bcell_program_celltype_summary.tsv")

fit_celltype_effect <- function(dat, cohort_name, cell_type_name) {
  work <- data.table::copy(dat[cohort == cohort_name & cell_type_analysis == cell_type_name])
  if (nrow(work) < 10L || length(unique(work$condition)) < 2L) return(NULL)
  work[, outcome_z := as.numeric(scale(mean_score))]
  work[, proportion_logit := qlogis(pmin(pmax(cell_proportion, 1e-4), 1 - 1e-4))]
  covariates <- available_covariates(work, c("age", "sex", "batch", "proportion_logit"))
  required <- c("outcome_z", "condition", covariates)
  work <- work[stats::complete.cases(work[, ..required])]
  if (nrow(work) < 10L || length(unique(work$condition)) < 2L) return(NULL)
  fit <- stats::lm(formula_from_terms("outcome_z", c("condition", covariates)), data = work)
  coef_table <- coef(summary(fit))
  if (!("conditionSLE" %in% rownames(coef_table))) return(NULL)
  estimate <- unname(coef_table["conditionSLE", "Estimate"])
  standard_error <- unname(coef_table["conditionSLE", "Std. Error"])
  data.table::data.table(cohort = cohort_name, cell_type_analysis = cell_type_name, estimate = estimate, standard_error = standard_error, conf_low = estimate - 1.96 * standard_error, conf_high = estimate + 1.96 * standard_error, statistic = estimate / standard_error, p_value = 2 * stats::pnorm(abs(estimate / standard_error), lower.tail = FALSE), n_profiles = nrow(work), n_donors = data.table::uniqueN(work$subject_id), covariates = paste(covariates, collapse = ";"))
}
bcell_cell_effects <- data.table::rbindlist(lapply(c("Child", "Adult"), function(cohort_name) {
  cell_types <- unique(bcell_scores[cohort == cohort_name, cell_type_analysis])
  data.table::rbindlist(lapply(cell_types, function(cell_type_name) fit_celltype_effect(bcell_scores, cohort_name, cell_type_name)), fill = TRUE)
}), fill = TRUE)
bcell_cell_effects[, fdr := stats::p.adjust(p_value, method = "BH"), by = cohort]
write_tsv(bcell_cell_effects, "01_bcell_localization/bcell_program_celltype_effects.tsv")

matching <- data.table::fread(input_paths[["program_matching"]], sep = "\t")
bcell_gene_list <- strsplit(matching[program == bcell_program, matched_gene_list][[1]], ";", fixed = TRUE)[[1]]
bcell_identity_genes <- intersect(bcell_gene_list, c("CD19", "CD22", "CD72", "CD79A", "CD79B", "MS4A1", "BANK1", "BLK", "BCL11A", "CR2", "CXCR5", "FCER2", "FCRL2", "FCRL3", "FCRL5", "GNG7", "HLA-DOA", "HLA-DOB", "POU2AF1", "SPIB", "TCL1A", "VPREB3", "AFF3", "BACH2", "BTLA"))
baff_nfk_genes <- intersect(bcell_gene_list, c("TNFSF13B", "BIRC3", "NFKB2", "RELB", "TNFRSF13C", "TNFRSF13B", "TNFRSF17"))
gene_sets <- list(full_B_cell_complement_program = bcell_gene_list, B_cell_identity = bcell_identity_genes, BAFF_NFkB_associated = baff_nfk_genes)
pseudobulk <- data.table::as.data.table(readRDS(input_paths[["pseudobulk"]]))
pseudobulk <- pseudobulk[comparison %in% c("Child_SLE_vs_HD", "Adult_SLE_vs_HD") & cell_label %in% c("B", "B_Naive")]
gene_decomposition <- data.table::rbindlist(lapply(names(gene_sets), function(set_name) {
  genes <- gene_sets[[set_name]]
  data.table::rbindlist(lapply(unique(pseudobulk$comparison), function(comparison_name) {
    data.table::rbindlist(lapply(unique(pseudobulk$cell_label), function(cell_label_name) {
      dat <- pseudobulk[comparison == comparison_name & cell_label == cell_label_name & gene %in% genes]
      if (!nrow(dat)) return(NULL)
      data.table::data.table(gene_set = set_name, comparison = comparison_name, cell_label = cell_label_name, n_genes_present = nrow(dat), mean_logFC = mean(dat$logFC, na.rm = TRUE), positive_fraction = mean(dat$logFC > 0, na.rm = TRUE), n_fdr_positive = sum(dat$FDR < 0.05, na.rm = TRUE), genes = paste(dat$gene, collapse = ";"))
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)
write_tsv(gene_decomposition, "01_bcell_localization/bcell_gene_set_decomposition.tsv")
gene_details <- pseudobulk[gene %in% bcell_gene_list, .(comparison, cell_label, gene, logFC, FDR)]
write_tsv(gene_details, "01_bcell_localization/bcell_program_gene_details.tsv")
 
