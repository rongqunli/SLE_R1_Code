if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
options(stringsAsFactors = FALSE, width = 120)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(sandwich)
  library(lmtest)
  library(irlba)
  library(ggforce)
  library(ggalluvial)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "analysis", "advanced_program_landscape.R"), winslash = "/", mustWork = TRUE)
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "advanced_program_landscape")
subdirs <- c(
  "00_audit", "01_donor_program_atlas", "02_multicellular_tensor",
  "03_trial_vectors", "figures", "logs", "checkpoints"
)
for (d in subdirs) dir.create(file.path(results_root, d), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(results_root, "logs", "analysis.log")
log_con <- file(log_file, open = "wt")
on.exit(close(log_con), add = TRUE)
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_con)
  flush(log_con)
}

write_tsv <- function(x, relative_path) {
  path <- file.path(results_root, relative_path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::as.data.table(x), path, sep = "\t", na = "NA", quote = FALSE)
  invisible(path)
}

atomic_save_rds <- function(x, path, compress = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(x, tmp, compress = compress)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not atomically save ", path)
  invisible(path)
}

save_figure <- function(plot, stem, width, height, dpi = 320) {
  png_path <- file.path(results_root, "figures", paste0(stem, ".png"))
  pdf_path <- file.path(results_root, "figures", paste0(stem, ".pdf"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE)
  invisible(c(png_path, pdf_path))
}

robust_lm_term <- function(formula, data, term) {
  fit <- stats::lm(formula, data = data)
  if (!term %in% names(stats::coef(fit))) stop("Missing model term: ", term)
  vc <- sandwich::vcovHC(fit, type = "HC3")
  ct <- lmtest::coeftest(fit, vcov. = vc)
  estimate <- unname(ct[term, "Estimate"])
  se <- unname(ct[term, "Std. Error"])
  p <- unname(ct[term, "Pr(>|t|)"])
  data.table(
    estimate = estimate, std_error = se,
    conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se,
    p_value = p, residual_df = stats::df.residual(fit), r_squared = summary(fit)$r.squared
  )
}

input_paths <- c(
  donor_scores = file.path(project_root, "pipeline", "3.5", "data", "processed", "donor_celltype_program_scores.tsv"),
  donor_proportions = file.path(project_root, "pipeline", "3.5", "data", "processed", "donor_celltype_proportions.tsv"),
  program_matching = file.path(project_root, "pipeline", "3.5", "data", "processed", "frozen_program_gene_matching.tsv"),
  pseudobulk_aggregations = "F:/program/Hei/SLE/GSE135779_pipeline/data/processed/checkpoints/05_pseudobulk/pseudobulk_aggregations.rds"
)
if (any(!file.exists(input_paths))) {
  stop("Missing inputs: ", paste(names(input_paths)[!file.exists(input_paths)], collapse = ", "))
}

input_manifest <- data.table(
  input_id = names(input_paths),
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  bytes = as.numeric(file.info(input_paths)$size),
  modified = as.character(file.info(input_paths)$mtime),
  md5 = unname(tools::md5sum(input_paths))
)
write_tsv(input_manifest, "00_audit/input_manifest.tsv")

package_names <- c(
  "data.table", "Matrix", "edgeR", "limma", "ggplot2", "patchwork",
  "sandwich", "lmtest", "irlba", "ggforce", "ggalluvial"
)
package_audit <- data.table(
  package = package_names,
  version = vapply(package_names, function(p) as.character(utils::packageVersion(p)), character(1)),
  implementation_role = c(
    "data orchestration", "sparse matrices", "TMM/logCPM normalization", "expression adjustment",
    "graphics", "figure composition", "HC3 covariance", "robust coefficient tests",
    "truncated SVD", "confidence ellipses", "multicellular flow diagram"
  )
)
package_audit <- rbind(
  package_audit,
  data.table(
    package = c("dreamlet", "scITD"), version = c(NA_character_, NA_character_),
    implementation_role = c(
      "not installed; donor-level HC3 program models and limma/edgeR pseudobulk normalization used",
      "not installed; auditable donor-mode HOSVD/Tucker core implemented with irlba"
    )
  ), fill = TRUE
)
write_tsv(package_audit, "00_audit/package_audit.tsv")

primary_programs <- c(
  "CYAN_NEUTROPHIL_IMMATURE", "GREENYELLOW_B_CELL_COMPLEMENT", "MAGENTA_IFN_ANTIVIRAL"
)
program_labels <- c(
  CYAN_NEUTROPHIL_IMMATURE = "Immature-neutrophil",
  GREENYELLOW_B_CELL_COMPLEMENT = "B-cell-enriched",
  MAGENTA_IFN_ANTIVIRAL = "IFN/antiviral"
)
program_colors <- c(
  "Immature-neutrophil" = "#D68A2D",
  "B-cell-enriched" = "#25837A",
  "IFN/antiviral" = "#C43C52"
)

theme_article <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = "#20242A"),
      plot.subtitle = element_text(color = "#555B65"),
      strip.background = element_rect(fill = "#F2F3F5", color = "#D0D3D8"),
      strip.text = element_text(face = "bold", color = "#252A31"),
      axis.text = element_text(color = "#30343B"),
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "#ECEEF1", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}




log_message("Starting donor-level all-cell-type program atlas")
scores <- fread(input_paths[["donor_scores"]], sep = "\t")
proportions <- fread(input_paths[["donor_proportions"]], sep = "\t")

scores <- scores[
  cohort == "Child" & analysis_set == "primary_ucell" &
    program %in% primary_programs & condition %in% c("HD", "SLE")
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

coverage <- scores[, .(
  n_donors = uniqueN(subject_id),
  n_samples = uniqueN(sample_id),
  median_cells = as.numeric(median(n_cells)),
  min_cells = as.numeric(min(n_cells))
), by = .(cell_type_analysis, condition)]
coverage_wide <- dcast(
  coverage, cell_type_analysis ~ condition,
  value.var = c("n_donors", "n_samples", "median_cells", "min_cells")
)
eligible_celltypes <- coverage_wide[
  n_donors_HD >= 10 & n_donors_SLE >= 20 & median_cells_HD >= 50 & median_cells_SLE >= 50,
  cell_type_analysis
]
coverage_wide[, eligible := cell_type_analysis %in% eligible_celltypes]
write_tsv(coverage_wide, "01_donor_program_atlas/celltype_coverage.tsv")

effect_results <- rbindlist(lapply(eligible_celltypes, function(cell_type_name) {
  rbindlist(lapply(primary_programs, function(program_name) {
    dat <- copy(scores[cell_type_analysis == cell_type_name & program == program_name])
    dat[, outcome_z := as.numeric(scale(mean_score))]
    covariates <- c("condition", "age_scaled")
    if (uniqueN(dat$sex) > 1) covariates <- c(covariates, "sex")
    if (uniqueN(dat$batch) > 1) covariates <- c(covariates, "batch")
    if (sd(dat$proportion_logit, na.rm = TRUE) > 0) covariates <- c(covariates, "proportion_logit")
    model_formula <- reformulate(covariates, response = "outcome_z")
    result <- robust_lm_term(model_formula, dat, "conditionSLE")
    result[, `:=`(
      cell_type = cell_type_name, program = program_name,
      n_HD = uniqueN(dat[condition == "HD", subject_id]),
      n_SLE = uniqueN(dat[condition == "SLE", subject_id]),
      model_formula = paste(deparse(model_formula), collapse = " "),
      effect_scale = "standardized donor mean UCell score; HC3 covariance"
    )]
    result
  }))
}))
effect_results[, fdr_global := p.adjust(p_value, method = "BH")]
effect_results[, fdr := p.adjust(p_value, method = "BH"), by = program]
effect_results[, significant := fdr < 0.05]
setcolorder(effect_results, c(
  "cell_type", "program", "estimate", "std_error", "conf_low", "conf_high",
  "p_value", "fdr", "fdr_global", "significant", "n_HD", "n_SLE", "model_formula",
  "effect_scale", "residual_df", "r_squared"
))
write_tsv(effect_results, "01_donor_program_atlas/donor_program_effects.tsv")

cell_order <- effect_results[program == "MAGENTA_IFN_ANTIVIRAL"][
  order(estimate), cell_type
]
atlas_plot_data <- copy(effect_results)
atlas_plot_data[, cell_type := factor(cell_type, levels = cell_order)]
atlas_plot_data[, program_label := factor(
  unname(program_labels[program]),
  levels = unname(program_labels[primary_programs])
)]
effect_limit <- max(abs(atlas_plot_data$estimate), na.rm = TRUE)

p_atlas <- ggplot(atlas_plot_data, aes(program_label, cell_type, fill = estimate)) +
  geom_tile(color = "white", linewidth = 1.1, width = 0.96, height = 0.94) +
  geom_tile(
    data = atlas_plot_data[significant == TRUE], fill = NA,
    color = "#111111", linewidth = 0.85, width = 0.96, height = 0.94
  ) +
  geom_text(aes(label = sprintf("%.2f", estimate)), size = 3.15, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "#F7F7F5", high = "#C53030", midpoint = 0,
    limits = c(-effect_limit, effect_limit), oob = squish,
    name = "Adjusted\nSLE effect"
  ) +
  coord_fixed(ratio = 0.72) +
  theme_article(10.5) +
  theme(
    axis.title = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), axis.text.x = element_text(angle = 24, hjust = 1),
    plot.margin = margin(8, 18, 8, 8)
  ) +
  labs(
    title = "Donor-level program effects across pediatric PBMC cell types",
    subtitle = "Composition-adjusted standardized SLE effects; black border: within-program FDR < 0.05"
  )
save_figure(p_atlas, "fig01_donor_program_effect_atlas", 9.8, 6.8)




log_message("Loading audited pseudobulk aggregations")
aggregations <- readRDS(input_paths[["pseudobulk_aggregations"]])
pb_counts <- aggregations$analysis$counts
pb_meta <- as.data.table(aggregations$analysis$metadata)
if (!identical(colnames(pb_counts), as.character(pb_meta$pseudobulk_id))) {
  stop("Pseudobulk count columns do not align to pseudobulk_id")
}

child_meta <- pb_meta[cohort == "Child" & condition %in% c("HD", "SLE")]
child_donors <- sort(unique(child_meta$subject_id))
candidate_tensor_celltypes <- c(
  "B_Naive", "CD14_Monocyte", "CD4_T", "CD4_T_Naive",
  "CD8_T_Cytotoxic", "CD8_T_GZMK", "CD8_T_Naive", "NK"
)
tensor_celltypes <- candidate_tensor_celltypes[
  vapply(candidate_tensor_celltypes, function(ct) {
    uniqueN(child_meta[cell_label == ct, subject_id]) == length(child_donors)
  }, logical(1))
]
if (length(child_donors) != 44L) stop("Expected 44 pediatric donors; observed ", length(child_donors))
if (length(tensor_celltypes) < 6L) stop("Insufficient complete cell types for tensor analysis")

tensor_meta <- child_meta[cell_label %in% tensor_celltypes]
tensor_meta[, subject_id := factor(subject_id, levels = child_donors)]
tensor_meta[, cell_label := factor(cell_label, levels = tensor_celltypes)]
setorder(tensor_meta, cell_label, subject_id)
tensor_counts <- pb_counts[, match(tensor_meta$pseudobulk_id, colnames(pb_counts)), drop = FALSE]
if (!identical(colnames(tensor_counts), as.character(tensor_meta$pseudobulk_id))) stop("Tensor subset alignment failed")

tensor_coverage <- tensor_meta[, .(
  donors = uniqueN(subject_id), median_cells = as.numeric(median(n_cells)), min_cells = as.numeric(min(n_cells)),
  HD = uniqueN(subject_id[condition == "HD"]), SLE = uniqueN(subject_id[condition == "SLE"])
), by = cell_label]
write_tsv(tensor_coverage, "02_multicellular_tensor/tensor_celltype_coverage.tsv")

program_matching <- fread(input_paths[["program_matching"]], sep = "\t")
program_gene_sets <- setNames(
  lapply(primary_programs, function(p) {
    strsplit(program_matching[program == p, matched_gene_list][[1]], ";", fixed = TRUE)[[1]]
  }),
  primary_programs
)

donor_meta <- unique(tensor_meta[, .(subject_id = as.character(subject_id), condition, age, sex, batch)])
setorder(donor_meta, subject_id)
if (!identical(donor_meta$subject_id, child_donors)) stop("Donor metadata ordering failed")
donor_meta[, condition := factor(condition, levels = c("HD", "SLE"))]
donor_meta[, age_scaled := as.numeric(scale(age))]

logcpm_list <- list()
adjusted_list <- list()
variance_list <- list()
for (ct in tensor_celltypes) {
  idx <- which(tensor_meta$cell_label == ct)
  ct_meta <- tensor_meta[idx]
  setorder(ct_meta, subject_id)
  ct_counts <- tensor_counts[, match(ct_meta$pseudobulk_id, colnames(tensor_counts)), drop = FALSE]
  if (!identical(as.character(ct_meta$subject_id), child_donors)) stop("Incomplete donor ordering for ", ct)
    dge <- edgeR::DGEList(counts = ct_counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
    colnames(logcpm) <- child_donors
  design_keep <- model.matrix(~ condition, donor_meta)
  covar <- model.matrix(~ age_scaled + sex, donor_meta)
  covar <- covar[, setdiff(colnames(covar), "(Intercept)"), drop = FALSE]
  adjusted <- limma::removeBatchEffect(
    logcpm, batch = donor_meta$batch,
    covariates = if (ncol(covar)) covar else NULL,
    design = design_keep
  )
  logcpm_list[[ct]] <- logcpm
  adjusted_list[[ct]] <- adjusted
  variance_list[[ct]] <- matrixStats::rowVars(adjusted)
}

variance_matrix <- do.call(cbind, variance_list)
colnames(variance_matrix) <- tensor_celltypes
detection_matrix <- sapply(tensor_celltypes, function(ct) {
  rowMeans(logcpm_list[[ct]] > 0)
})
eligible_genes <- rownames(variance_matrix)[rowSums(detection_matrix >= 0.20) >= ceiling(length(tensor_celltypes) / 2)]
median_variance <- apply(variance_matrix[eligible_genes, , drop = FALSE], 1, median, na.rm = TRUE)
selected_genes <- names(sort(median_variance, decreasing = TRUE))[seq_len(min(2500L, length(median_variance)))]

feature_blocks <- list()
feature_table <- list()
for (ct in tensor_celltypes) {
  mat <- adjusted_list[[ct]][selected_genes, child_donors, drop = FALSE]
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  feature_blocks[[ct]] <- t(z)
  feature_table[[ct]] <- data.table(
    feature_index_within_celltype = seq_along(selected_genes),
    cell_type = ct, gene = selected_genes
  )
}
unfold <- do.call(cbind, feature_blocks)
features <- rbindlist(feature_table)
if (nrow(unfold) != 44L || ncol(unfold) != nrow(features)) stop("Tensor unfolding dimension mismatch")

n_factors <- min(10L, nrow(unfold) - 2L)
set.seed(20260809)
sv <- irlba::irlba(unfold, nv = n_factors, nu = n_factors, center = FALSE, scale = FALSE)
donor_factor_scores <- sweep(sv$u, 2, sv$d, `*`)
feature_loadings <- sv$v
factor_names <- paste0("Factor", seq_len(n_factors))
colnames(donor_factor_scores) <- factor_names
colnames(feature_loadings) <- factor_names
rownames(donor_factor_scores) <- child_donors

global_program_scores <- scores[
  cell_type_analysis %in% tensor_celltypes,
  .(global_score = mean(mean_score, na.rm = TRUE)),
  by = .(subject_id, program)
]
global_program_wide <- dcast(global_program_scores, subject_id ~ program, value.var = "global_score")
setorder(global_program_wide, subject_id)
if (!identical(global_program_wide$subject_id, child_donors)) stop("Global program score donor alignment failed")

factor_program_alignment <- rbindlist(lapply(seq_len(n_factors), function(k) {
  rbindlist(lapply(primary_programs, function(p) {
    rho <- suppressWarnings(cor(donor_factor_scores[, k], global_program_wide[[p]], method = "spearman"))
    data.table(factor = factor_names[[k]], program = p, spearman_rho = rho)
  }))
}))
ifn_factor <- factor_program_alignment[program == "MAGENTA_IFN_ANTIVIRAL"][
  which.max(abs(spearman_rho)), factor
]
ifn_factor_index <- match(ifn_factor, factor_names)
ifn_orientation <- sign(factor_program_alignment[
  factor == ifn_factor & program == "MAGENTA_IFN_ANTIVIRAL", spearman_rho
])
if (!is.finite(ifn_orientation) || ifn_orientation == 0) ifn_orientation <- 1
donor_factor_scores[, ifn_factor_index] <- donor_factor_scores[, ifn_factor_index] * ifn_orientation
feature_loadings[, ifn_factor_index] <- feature_loadings[, ifn_factor_index] * ifn_orientation
factor_program_alignment[factor == ifn_factor, spearman_rho := spearman_rho * ifn_orientation]
donor_factor_scores <- scale(donor_factor_scores)
colnames(donor_factor_scores) <- factor_names
rownames(donor_factor_scores) <- child_donors

factor_associations <- rbindlist(lapply(seq_len(n_factors), function(k) {
  dat <- copy(donor_meta)
  dat[, factor_score := donor_factor_scores[, k]]
  result <- robust_lm_term(
    factor_score ~ condition + age_scaled + sex + batch,
    dat, "conditionSLE"
  )
  result[, factor := factor_names[[k]]]
  result
}))
factor_associations[, fdr := p.adjust(p_value, method = "BH")]
factor_associations[, explained_fraction := sv$d^2 / sum(unfold^2)]
factor_associations[, selected_ifn_aligned := factor == ifn_factor]
write_tsv(factor_associations, "02_multicellular_tensor/factor_disease_associations.tsv")
write_tsv(factor_program_alignment, "02_multicellular_tensor/factor_program_alignment.tsv")

loading_long <- cbind(features, as.data.table(feature_loadings))
write_tsv(
  loading_long[, c("cell_type", "gene", factor_names), with = FALSE],
  "02_multicellular_tensor/all_factor_loadings.tsv"
)

ifn_loading <- loading_long[, .(
  cell_type, gene, loading = get(ifn_factor),
  is_ifn_gene = gene %in% program_gene_sets[["MAGENTA_IFN_ANTIVIRAL"]],
  is_bcell_gene = gene %in% program_gene_sets[["GREENYELLOW_B_CELL_COMPLEMENT"]],
  is_neutrophil_gene = gene %in% program_gene_sets[["CYAN_NEUTROPHIL_IMMATURE"]]
)]
write_tsv(ifn_loading[order(-loading)], "02_multicellular_tensor/selected_ifn_factor_loadings.tsv")

gene_set_loading_summary <- rbindlist(lapply(primary_programs, function(p) {
  genes <- program_gene_sets[[p]]
  ifn_loading[, .(
    program = p,
    detected_features = sum(gene %in% genes),
    mean_loading = mean(loading[gene %in% genes], na.rm = TRUE),
    mean_abs_loading = mean(abs(loading[gene %in% genes]), na.rm = TRUE),
    background_mean_abs_loading = mean(abs(loading[!gene %in% genes]), na.rm = TRUE),
    wilcox_p = tryCatch(
      wilcox.test(abs(loading[gene %in% genes]), abs(loading[!gene %in% genes]))$p.value,
      error = function(e) NA_real_
    )
  )]
}))
gene_set_loading_summary[, fdr := p.adjust(wilcox_p, method = "BH")]
write_tsv(gene_set_loading_summary, "02_multicellular_tensor/selected_factor_gene_set_enrichment.tsv")

donor_scores_out <- cbind(
  donor_meta,
  as.data.table(donor_factor_scores)
)
write_tsv(donor_scores_out, "02_multicellular_tensor/donor_factor_scores.tsv")


log_message("Running tensor loading bootstrap stability")
set.seed(20260810)
n_boot <- 100L
full_loading <- feature_loadings[, ifn_factor_index]
bootstrap_stability <- rbindlist(lapply(seq_len(n_boot), function(b) {
  idx <- sample(seq_len(nrow(unfold)), replace = TRUE)
  boot_sv <- irlba::irlba(unfold[idx, , drop = FALSE], nv = n_factors, nu = n_factors)
  cors <- vapply(seq_len(n_factors), function(k) {
    suppressWarnings(cor(full_loading, boot_sv$v[, k]))
  }, numeric(1))
  best <- which.max(abs(cors))
  data.table(
    bootstrap = b, matched_factor = best,
    loading_correlation = abs(cors[[best]]),
    signed_correlation = cors[[best]]
  )
}))
write_tsv(bootstrap_stability, "02_multicellular_tensor/tensor_loading_bootstrap_stability.tsv")

selected_factor_data <- data.table(
  subject_id = child_donors,
  condition = donor_meta$condition,
  score = donor_factor_scores[, ifn_factor_index]
)
p_factor_score <- ggplot(selected_factor_data, aes(condition, score, color = condition, fill = condition)) +
  geom_violin(width = 0.78, trim = TRUE, alpha = 0.14, linewidth = 0) +
  geom_boxplot(width = 0.23, outlier.shape = NA, fill = "white", color = "#34383E", linewidth = 0.55) +
  geom_point(position = position_jitter(width = 0.09, seed = 20260809), size = 2.1, alpha = 0.9) +
  scale_color_manual(values = c(HD = "#2C6E9B", SLE = "#C23B3B")) +
  scale_fill_manual(values = c(HD = "#2C6E9B", SLE = "#C23B3B")) +
  theme_article(10.5) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  labs(
    title = paste0(ifn_factor, ": donor scores"),
    subtitle = "Factor orientation fixed to positive IFN-program alignment",
    y = "Standardized multicellular factor score"
  )

alignment_plot_data <- copy(factor_program_alignment)
alignment_plot_data[, factor := factor(factor, levels = rev(factor_names))]
alignment_plot_data[, program_label := factor(
  unname(program_labels[program]), levels = unname(program_labels[primary_programs])
)]
p_alignment <- ggplot(alignment_plot_data, aes(program_label, factor, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2.8) +
  geom_rect(
    data = alignment_plot_data[factor == ifn_factor],
    aes(xmin = as.numeric(program_label) - 0.48, xmax = as.numeric(program_label) + 0.48,
        ymin = as.numeric(factor) - 0.48, ymax = as.numeric(factor) + 0.48),
    inherit.aes = FALSE, fill = NA, color = "#111111", linewidth = 0.7
  ) +
  scale_fill_gradient2(low = "#2B6CB0", mid = "white", high = "#C53030", midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
  theme_article(9.5) +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 24, hjust = 1), panel.grid = element_blank()) +
  labs(title = "Factor alignment with frozen donor programs")

gene_rank <- ifn_loading[, .(mean_loading = mean(loading)), by = gene][order(-mean_loading)]
heatmap_genes <- unique(c(
  gene_rank[gene %in% program_gene_sets[["MAGENTA_IFN_ANTIVIRAL"]], head(gene, 18)],
  intersect(c("STAT1", "STAT2", "IRF7", "JAK2", "SOCS1"), selected_genes)
))
heatmap_data <- ifn_loading[gene %in% heatmap_genes]
gene_levels <- gene_rank[gene %in% heatmap_genes, gene]
heatmap_data[, gene := factor(gene, levels = rev(gene_levels))]
heatmap_data[, cell_type := factor(cell_type, levels = tensor_celltypes)]
heat_limit <- stats::quantile(abs(heatmap_data$loading), 0.98, na.rm = TRUE)
p_loading <- ggplot(heatmap_data, aes(cell_type, gene, fill = loading)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "#F7F7F5", high = "#C53030", midpoint = 0,
    limits = c(-heat_limit, heat_limit), oob = squish, name = "Factor\nloading"
  ) +
  theme_article(9.3) +
  theme(
    axis.title = element_blank(), axis.text.x = element_text(angle = 32, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "IFN-aligned factor loadings across cell types")

p_stability <- ggplot(bootstrap_stability, aes(loading_correlation)) +
  geom_histogram(bins = 18, fill = "#7566A8", color = "white", linewidth = 0.45) +
  geom_vline(
    xintercept = median(bootstrap_stability$loading_correlation),
    color = "#252A31", linewidth = 0.8, linetype = 2
  ) +
  annotate(
    "label", x = median(bootstrap_stability$loading_correlation), y = Inf,
    label = sprintf("median = %.2f", median(bootstrap_stability$loading_correlation)),
    vjust = 1.3, hjust = 1.05, size = 3.1, fill = "white"
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_article(10) +
  labs(
    title = "Bootstrap loading stability",
    subtitle = paste0(n_boot, " donor-resampled decompositions"),
    x = "Maximum absolute loading correlation", y = "Bootstrap runs"
  )

tensor_combined <- (p_factor_score | p_alignment) / (p_loading | p_stability) +
  plot_annotation(
    title = "Donor-resolved multicellular factorization of the pediatric PBMC atlas",
    subtitle = paste0(
      "44 donors × ", length(tensor_celltypes), " complete cell types × ",
      length(selected_genes), " variable genes; donor labels were not used during factorization"
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 10.5, color = "#555B65"))
  )
save_figure(tensor_combined, "fig02_multicellular_tensor_factor", 15.5, 10.8)

flow_data <- ifn_loading[loading > 0][order(cell_type, -loading), head(.SD, 3), by = cell_type]
flow_data[, target := paste0(ifn_factor, "\nIFN-aligned")]
flow_data[, weight := abs(loading)]
flow_data[, cell_type := factor(cell_type, levels = tensor_celltypes)]
p_flow <- ggplot(
  flow_data,
  aes(axis1 = cell_type, axis2 = gene, axis3 = target, y = weight, fill = cell_type)
) +
  geom_alluvium(alpha = 0.72, width = 0.18, knot.pos = 0.42) +
  geom_stratum(width = 0.18, fill = "#F4F5F6", color = "#50555C", linewidth = 0.35) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.8, color = "#20242A") +
  scale_x_discrete(limits = c("Cell type", "Top loading genes", "Multicellular factor"), expand = c(0.08, 0.08)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "#555B65")) +
  labs(
    title = "Cell-type contributions to the IFN-aligned multicellular factor",
    subtitle = "Flows display the three highest positive gene loadings per complete cell type; widths encode absolute loading"
  )
save_figure(p_flow, "fig03_multicellular_factor_flow", 13.5, 7.2)

atomic_save_rds(
  list(
    donor_factor_scores = donor_scores_out,
    factor_program_alignment = factor_program_alignment,
    factor_associations = factor_associations,
    selected_ifn_factor = ifn_factor,
    selected_factor_loadings = ifn_loading,
    singular_values = sv$d,
    selected_genes = selected_genes,
    cell_types = tensor_celltypes
  ),
  file.path(results_root, "checkpoints", "multicellular_tensor_results.rds"),
  compress = TRUE
)
rm(aggregations, pb_counts, tensor_counts, logcpm_list, adjusted_list, feature_blocks, unfold)
gc()
