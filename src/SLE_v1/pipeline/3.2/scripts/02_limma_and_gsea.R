if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
project_root <- Sys.getenv("SLE_PROJECT_ROOT")
stopifnot(nzchar(project_root))
analysis_root <- file.path(project_root, "pipeline", "3.2")
results_root <- file.path(project_root, "results", "3.2")

required <- c("data.table", "limma", "fgsea", "ggplot2", "ComplexHeatmap", "circlize")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

cfg_tab <- data.table::fread(file.path(analysis_root, "config", "analysis_contract.tsv"))
cfg <- setNames(as.list(cfg_tab[["value"]]), cfg_tab[["key"]])
set.seed(as.integer(cfg$random_seed))
resolve_resource <- function(path) {
  if (grepl("^[A-Za-z]:[/\\\\]", path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  normalizePath(file.path(analysis_root, path), winslash = "/", mustWork = TRUE)
}

object <- readRDS(file.path(results_root, "01_expression", "full_gene_expression.rds"))
expr <- object$expression
meta <- data.table::as.data.table(object$metadata)
meta[, condition := factor(condition, levels = c(cfg$control_label, cfg$case_label))]

out_dir <- file.path(results_root, "02_deg_gsea")
fig_dir <- file.path(results_root, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

design <- stats::model.matrix(~ condition, data = meta)
fit <- limma::eBayes(limma::lmFit(expr, design), robust = TRUE, trend = TRUE)
coef_name <- grep("^condition", colnames(design), value = TRUE)
if (length(coef_name) != 1L) stop("Could not identify the case-control coefficient.")
deg <- data.table::as.data.table(limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none"), keep.rownames = "gene")
data.table::setnames(deg, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"), c("logFC", "mean_expression", "t_stat", "p_value", "fdr", "lods"))
deg[, direction := data.table::fcase(
  fdr <= as.numeric(cfg$primary_fdr) & logFC >= as.numeric(cfg$primary_abs_logfc), "SLE_up",
  fdr <= as.numeric(cfg$primary_fdr) & logFC <= -as.numeric(cfg$primary_abs_logfc), "SLE_down",
  default = "not_primary"
)]

balanced_once <- function(i) {
  control <- which(meta$condition == cfg$control_label)
  cases <- sample(which(meta$condition == cfg$case_label), length(control), replace = FALSE)
  keep <- c(control, cases)
  d <- droplevels(meta$condition[keep])
  x <- stats::model.matrix(~ d)
  f <- limma::eBayes(limma::lmFit(expr[, keep, drop = FALSE], x), robust = TRUE, trend = TRUE)
  data.table::data.table(gene = rownames(expr), iteration = i, logFC = f$coefficients[, 2], p_value = f$p.value[, 2])
}
balanced <- data.table::rbindlist(lapply(seq_len(as.integer(cfg$balanced_resamples)), balanced_once))
stability <- balanced[, .(
  balanced_median_logFC = stats::median(logFC, na.rm = TRUE),
  direction_agreement = max(mean(logFC > 0, na.rm = TRUE), mean(logFC < 0, na.rm = TRUE)),
  nominal_hit_fraction = mean(p_value < 0.05, na.rm = TRUE)
), by = gene]
deg <- merge(deg, stability, by = "gene", all.x = TRUE)
data.table::fwrite(deg, file.path(out_dir, "full_transcriptome_limma.tsv"), sep = "\t")

run_sensitivity_model <- function(model_name, keep, covariates = character()) {
  sub_meta <- droplevels(as.data.frame(meta[keep]))
  rhs <- paste(c("condition", covariates), collapse = " + ")
  design_s <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data = sub_meta)
  if (qr(design_s)$rank < ncol(design_s)) stop("Rank-deficient sensitivity design: ", model_name)
  coef_s <- grep("^condition", colnames(design_s), value = TRUE)
  if (length(coef_s) != 1L) stop("Condition coefficient not identifiable in: ", model_name)
  fit_s <- limma::eBayes(limma::lmFit(expr[, keep, drop = FALSE], design_s), robust = TRUE, trend = TRUE)
  tab <- data.table::as.data.table(limma::topTable(fit_s, coef = coef_s, number = Inf, sort.by = "none"), keep.rownames = "gene")
  tab[, .(
    model = model_name, gene, logFC, t_stat = t, p_value = P.Value, fdr = adj.P.Val,
    n_samples = sum(keep), n_sle = sum(sub_meta$condition == cfg$case_label),
    n_healthy = sum(sub_meta$condition == cfg$control_label)
  )]
}

trial_batches <- c("ILLUMINATE-1", "ILLUMINATE-2")
sensitivity_models <- data.table::rbindlist(list(
  run_sensitivity_model("trial_batch_adjusted", meta$batch %in% trial_batches, "batch"),
  run_sensitivity_model("female_only", meta$sex == "Female" & !is.na(meta$sex))
))
primary_compare <- deg[, .(gene, primary_logFC = logFC, primary_fdr = fdr, primary_direction = direction)]
sensitivity_models <- merge(sensitivity_models, primary_compare, by = "gene", all.x = TRUE)
sensitivity_models[, direction_matches_primary := sign(logFC) == sign(primary_logFC)]
sensitivity_summary <- sensitivity_models[, .(
  n_samples = unique(n_samples), n_sle = unique(n_sle), n_healthy = unique(n_healthy),
  spearman_logFC = stats::cor(logFC, primary_logFC, method = "spearman", use = "complete.obs"),
  all_gene_direction_agreement = mean(direction_matches_primary, na.rm = TRUE),
  primary_gene_direction_agreement = mean(direction_matches_primary[primary_direction != "not_primary"], na.rm = TRUE),
  primary_gene_fdr_replication = mean(fdr[primary_direction != "not_primary"] <= as.numeric(cfg$primary_fdr), na.rm = TRUE)
), by = model]
covariate_audit <- data.table::data.table(
  covariate = c("age", "sex", "race", "batch"),
  healthy_available = c(
    sum(meta$condition == cfg$control_label & is.finite(suppressWarnings(as.numeric(meta$age)))),
    sum(meta$condition == cfg$control_label & !is.na(meta$sex) & nzchar(meta$sex)),
    sum(meta$condition == cfg$control_label & !is.na(meta$race) & nzchar(meta$race)),
    sum(meta$condition == cfg$control_label & !is.na(meta$batch) & nzchar(meta$batch))
  ),
  healthy_total = sum(meta$condition == cfg$control_label),
  decision = c(
    "not_estimable_no_healthy_values",
    "female_restricted_sensitivity",
    "not_estimable_no_healthy_values",
    "trial_stratified_adjustment"
  )
)
data.table::fwrite(sensitivity_models, file.path(out_dir, "covariate_sensitivity_gene_results.tsv"), sep = "\t")
data.table::fwrite(sensitivity_summary, file.path(out_dir, "covariate_sensitivity_summary.tsv"), sep = "\t")
data.table::fwrite(covariate_audit, file.path(out_dir, "covariate_estimability_audit.tsv"), sep = "\t")

read_gmt <- function(path, prefix) {
  if (!file.exists(path)) return(list())
  x <- fgsea::gmtPathways(path)
  names(x) <- paste0(prefix, "::", names(x))
  x
}
pathways <- c(
  read_gmt(resolve_resource(cfg$go_gmt), "GO"),
  read_gmt(resolve_resource(cfg$kegg_gmt), "KEGG")
)

program_file <- resolve_resource(cfg$program_definitions)
programs <- data.table::fread(program_file)
if (all(c("program", "genes") %in% names(programs))) {
  custom_long <- programs[, .(gene = unlist(strsplit(genes, ";", fixed = TRUE))), by = program]
  custom <- split(toupper(trimws(custom_long$gene)), paste0("FROZEN::", custom_long$program))
  pathways <- c(pathways, custom)
} else if (all(c("program", "gene") %in% names(programs))) {
  custom <- split(toupper(trimws(programs$gene)), paste0("FROZEN::", programs$program))
  pathways <- c(pathways, custom)
}
pathways <- lapply(pathways, unique)
if (!length(pathways)) stop("No gene sets were available for GSEA.")

ranks <- deg$t_stat
names(ranks) <- deg$gene
ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
gsea <- fgsea::fgseaMultilevel(
  pathways = pathways, stats = ranks,
  minSize = as.integer(cfg$min_gene_set_size),
  maxSize = as.integer(cfg$max_gene_set_size),
  eps = as.numeric(cfg$fgsea_eps)
)
gsea <- data.table::as.data.table(gsea)
gsea[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
gsea[, collection := sub("::.*$", "", pathway)]
gsea[, label := sub("^[^:]+::", "", pathway)]
gsea[, abs_NES := abs(NES)]
data.table::setorder(gsea, padj, -abs_NES)
data.table::fwrite(gsea, file.path(out_dir, "fgsea_all_collections.tsv"), sep = "\t")

gsea_sensitivity <- data.table::rbindlist(lapply(unique(sensitivity_models$model), function(model_name) {
  model_stats <- sensitivity_models[model == model_name]
  ranks_s <- model_stats$t_stat
  names(ranks_s) <- model_stats$gene
  ranks_s <- sort(ranks_s[is.finite(ranks_s) & !duplicated(names(ranks_s))], decreasing = TRUE)
  ans <- fgsea::fgseaMultilevel(
    pathways = pathways, stats = ranks_s,
    minSize = as.integer(cfg$min_gene_set_size),
    maxSize = as.integer(cfg$max_gene_set_size),
    eps = as.numeric(cfg$fgsea_eps)
  )
  ans <- data.table::as.data.table(ans)
  ans[, `:=`(
    model = model_name,
    leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"),
    collection = sub("::.*$", "", pathway),
    label = sub("^[^:]+::", "", pathway),
    abs_NES = abs(NES)
  )]
  ans
}), use.names = TRUE, fill = TRUE)
data.table::setorder(gsea_sensitivity, model, padj, -abs_NES)
primary_pathway <- gsea[, .(pathway, primary_NES = NES, primary_padj = padj)]
gsea_sensitivity <- merge(gsea_sensitivity, primary_pathway, by = "pathway", all.x = TRUE)
gsea_sensitivity[, direction_matches_primary := sign(NES) == sign(primary_NES)]
gsea_sensitivity_summary <- gsea_sensitivity[, .(
  tested_pathways = .N,
  direction_agreement = mean(direction_matches_primary, na.rm = TRUE),
  primary_significant_pathways = sum(primary_padj <= as.numeric(cfg$primary_fdr), na.rm = TRUE),
  primary_significant_direction_agreement = mean(direction_matches_primary[primary_padj <= as.numeric(cfg$primary_fdr)], na.rm = TRUE),
  primary_significant_fdr_replication = mean(padj[primary_padj <= as.numeric(cfg$primary_fdr)] <= as.numeric(cfg$primary_fdr), na.rm = TRUE)
), by = .(model, collection)]
data.table::fwrite(gsea_sensitivity, file.path(out_dir, "fgsea_covariate_sensitivity.tsv"), sep = "\t")
data.table::fwrite(gsea_sensitivity_summary, file.path(out_dir, "fgsea_covariate_sensitivity_summary.tsv"), sep = "\t")

volcano <- ggplot2::ggplot(deg, ggplot2::aes(logFC, -log10(pmax(fdr, 1e-300)), color = direction)) +
  ggplot2::geom_point(alpha = 0.38, size = 0.8) +
  ggplot2::scale_color_manual(values = c(SLE_up = "#C73E1D", SLE_down = "#2166AC", not_primary = "#B8B8B8")) +
  ggplot2::geom_vline(xintercept = c(-1, 1) * as.numeric(cfg$primary_abs_logfc), linetype = 3, color = "#666666") +
  ggplot2::labs(x = "log2 fold change (SLE vs healthy)", y = "-log10 FDR", color = NULL, title = "GSE88884 full-transcriptome differential landscape") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top")
if (requireNamespace("ggrepel", quietly = TRUE)) {
  label_data <- data.table::rbindlist(list(
    deg[direction == "SLE_up"][order(-logFC)][1:10],
    deg[direction == "SLE_down"][order(logFC)][1:6]
  ))
  volcano <- volcano + ggrepel::geom_text_repel(
    data = label_data, ggplot2::aes(label = gene),
    size = 3.1, color = "#222222", box.padding = 0.35, point.padding = 0.2,
    min.segment.length = 0, max.overlaps = Inf, seed = as.integer(cfg$random_seed)
  )
}
ggplot2::ggsave(file.path(fig_dir, "fig01A_full_transcriptome_volcano.png"), volcano, width = 8.4, height = 6.2, dpi = 320)

if (nrow(gsea)) {
  top_gsea <- gsea[is.finite(padj) & (NES > 0 | collection == "FROZEN")][
    order(padj, -NES), head(.SD, if (collection[[1]] == "FROZEN") 3L else 14L), by = collection
  ]
  top_gsea[, label := factor(label, levels = rev(unique(label[order(NES)])))]
  gsea_dot <- ggplot2::ggplot(top_gsea, ggplot2::aes(NES, label, size = -log10(pmax(padj, 1e-12)), color = NES)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
    ggplot2::labs(x = "Normalized enrichment score", y = NULL, size = "-log10 FDR", color = "NES", title = "Pathway enrichment landscape") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(file.path(fig_dir, "fig01B_gsea_landscape.png"), gsea_dot, width = 12.5, height = 10, dpi = 320)
} else {
  warning("GSEA returned no eligible pathways; enrichment figures were skipped.")
}

sig <- gsea[padj <= 0.05][order(padj, -abs(NES))]
if (nrow(sig)) {
  chosen <- unique(data.table::rbindlist(list(
    sig[collection == "FROZEN" & NES > 0],
    sig[collection == "GO" & NES > 0][1:3],
    sig[collection == "KEGG" & NES > 0][1:3]
  )), by = "pathway")
  running <- lapply(chosen$pathway, function(p) {
    pdat <- fgsea::plotEnrichment(pathways[[p]], ranks) +
      ggplot2::labs(title = sub("^[^:]+::", "", p), subtitle = p) +
      ggplot2::theme_minimal(base_size = 9)
    pdat
  })
  if (length(running) && requireNamespace("patchwork", quietly = TRUE)) {
    panel <- patchwork::wrap_plots(running, ncol = 2)
    ggplot2::ggsave(file.path(fig_dir, "fig01C_gsea_running_scores.png"), panel, width = 13, height = 3.2 * ceiling(length(running) / 2), dpi = 300, limitsize = FALSE)
  }
}

saveRDS(
  list(
    deg = deg, gsea = gsea, gsea_sensitivity = gsea_sensitivity,
    sensitivity_models = sensitivity_models, pathways = pathways, ranks = ranks
  ),
  file.path(out_dir, "deg_gsea_results.rds")
)
