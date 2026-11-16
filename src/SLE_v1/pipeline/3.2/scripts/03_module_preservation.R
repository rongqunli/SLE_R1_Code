if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
project_root <- Sys.getenv("SLE_PROJECT_ROOT")
stopifnot(nzchar(project_root))
analysis_root <- file.path(project_root, "pipeline", "3.2")
results_root <- file.path(project_root, "results", "3.2")
required <- c("data.table", "WGCNA", "ggplot2", "doParallel", "foreach")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

cfg_tab <- data.table::fread(file.path(analysis_root, "config", "analysis_contract.tsv"))
cfg <- setNames(as.list(cfg_tab[["value"]]), cfg_tab[["key"]])
resolve <- function(x) normalizePath(file.path(analysis_root, x), winslash = "/", mustWork = TRUE)
set.seed(as.integer(cfg$random_seed))
WGCNA::allowWGCNAThreads()

ref_expr_obj <- readRDS(resolve(cfg$reference_expression))
ref_wgcna <- readRDS(resolve(cfg$reference_wgcna))
test_obj <- readRDS(resolve(cfg$test_expression))

ref_expr <- ref_expr_obj$expression
test_expr <- test_obj$expression
ref_modules <- data.table::as.data.table(ref_wgcna$module_assignments)
if (!all(c("gene", "module") %in% names(ref_modules))) stop("Reference module assignment contract failed.")
ref_modules[, gene := toupper(gene)]
rownames(ref_expr) <- toupper(rownames(ref_expr))
rownames(test_expr) <- toupper(rownames(test_expr))

ref_samples <- intersect(ref_wgcna$analysis_data$sample_id, colnames(ref_expr))
test_meta <- data.table::as.data.table(test_obj$metadata)
test_samples <- test_meta[condition == "SLE", sample_id]
test_samples <- intersect(test_samples, colnames(test_expr))
common <- Reduce(intersect, list(ref_modules$gene, rownames(ref_expr), rownames(test_expr)))
if (length(common) < 1000L) stop("Fewer than 1,000 common module genes are available across cohorts.")

module_map <- ref_modules[match(common, gene)]
keep <- !is.na(module_map$module) & module_map$module != "grey"
common <- common[keep]
module_colors <- module_map$module[keep]

ref_dat <- t(ref_expr[common, ref_samples, drop = FALSE])
test_dat <- t(test_expr[common, test_samples, drop = FALSE])
good_ref <- WGCNA::goodSamplesGenes(ref_dat, verbose = 0)
good_test <- WGCNA::goodSamplesGenes(test_dat, verbose = 0)
good_gene <- good_ref$goodGenes & good_test$goodGenes
ref_dat <- ref_dat[good_ref$goodSamples, good_gene, drop = FALSE]
test_dat <- test_dat[good_test$goodSamples, good_gene, drop = FALSE]
module_colors <- module_colors[good_gene]
common <- common[good_gene]
colnames(ref_dat) <- colnames(test_dat) <- common

multi_expr <- list()
multi_expr[[cfg$reference_name]] <- list(data = ref_dat)
multi_expr[[cfg$test_name]] <- list(data = test_dat)
color_list <- list()
color_list[[cfg$reference_name]] <- module_colors

worker_count <- max(1L, min(as.integer(cfg$parallel_workers), parallel::detectCores(logical = FALSE)))
cluster_backend <- parallel::makeCluster(worker_count)
doParallel::registerDoParallel(cluster_backend)
pres <- tryCatch(
  WGCNA::modulePreservation(
    multiData = multi_expr,
    multiColor = color_list,
    referenceNetworks = 1,
    nPermutations = as.integer(cfg$n_permutations),
    randomSeed = as.integer(cfg$random_seed),
    networkType = cfg$network_type,
    corFnc = cfg$cor_function,
    verbose = 3,
    quickCor = as.numeric(cfg$quick_cor),
    parallelCalculation = TRUE,
    savePermutedStatistics = FALSE
  ),
  finally = {
    parallel::stopCluster(cluster_backend)
    foreach::registerDoSEQ()
  }
)
z <- data.table::as.data.table(pres$preservation$Z[[1]][[2]], keep.rownames = "module")
obs <- data.table::as.data.table(pres$preservation$observed[[1]][[2]], keep.rownames = "module")
preservation <- merge(z, obs, by = "module", all = TRUE, suffixes = c(".z", ".observed"))
z_col <- grep("Zsummary", names(preservation), value = TRUE)[1]
rank_col <- grep("medianRank", names(preservation), value = TRUE)[1]
if (is.na(z_col) || is.na(rank_col)) stop("WGCNA preservation output did not expose Zsummary and medianRank.")
preservation[, preservation_class := data.table::fcase(
  get(z_col) >= as.numeric(cfg$zsummary_strong), "strong",
  get(z_col) >= as.numeric(cfg$zsummary_moderate), "moderate",
  default = "weak_or_none"
)]

test_me_obj <- WGCNA::moduleEigengenes(test_dat, colors = module_colors, softPower = 1, excludeGrey = TRUE)
test_mes <- WGCNA::orderMEs(test_me_obj$eigengenes)
test_kme <- WGCNA::signedKME(test_dat, test_mes, outputColumnName = "kME")
test_membership <- data.table::data.table(gene = common, module = module_colors)
for (m in unique(module_colors)) {
  col <- paste0("kME", m)
  test_membership[module == m, adult_kME := if (col %in% colnames(test_kme)) test_kme[gene, col] else NA_real_]
}
ref_membership <- data.table::as.data.table(ref_wgcna$membership)[, .(gene = toupper(gene), module, child_kME = kME, child_hub_rank = hub_rank)]
membership <- merge(ref_membership, test_membership, by = c("gene", "module"), all.x = TRUE)
membership[, adult_hub_rank := data.table::frank(-abs(adult_kME), ties.method = "average"), by = module]
membership_concordance <- membership[is.finite(child_kME) & is.finite(adult_kME), .(
  n_genes = .N,
  spearman_kME = stats::cor(child_kME, adult_kME, method = "spearman"),
  spearman_abs_kME = stats::cor(abs(child_kME), abs(adult_kME), method = "spearman")
), by = module]

test_me <- data.table::as.data.table(test_mes, keep.rownames = "sample_id")
test_meta_small <- unique(test_meta[, .(sample_id, SLEDAI, C3, C4, anti_dsDNA, age, sex, race, batch)])
clinical <- merge(test_me, test_meta_small, by = "sample_id", all.x = TRUE)
traits <- c("SLEDAI", "C3", "C4", "anti_dsDNA")
me_cols <- grep("^ME", names(clinical), value = TRUE)
clinical_assoc <- data.table::rbindlist(lapply(me_cols, function(me) {
  data.table::rbindlist(lapply(traits, function(trait) {
    dat <- clinical[, .(y = suppressWarnings(as.numeric(get(me))), x = suppressWarnings(as.numeric(get(trait))))]
    dat <- dat[is.finite(x) & is.finite(y)]
    if (nrow(dat) < 30L || stats::sd(dat$x) == 0) return(NULL)
    fit <- stats::lm(scale(y) ~ scale(x), data = dat)
    co <- summary(fit)$coefficients[2, ]
    data.table::data.table(module = sub("^ME", "", me), trait = trait, n = nrow(dat), beta = co[1], se = co[2], p_value = co[4])
  }), fill = TRUE)
}), fill = TRUE)
if (nrow(clinical_assoc)) clinical_assoc[, fdr := stats::p.adjust(p_value, method = "BH")]

out_dir <- file.path(results_root, "03_preservation")
fig_dir <- file.path(results_root, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(preservation, file.path(out_dir, "module_preservation_statistics.tsv"), sep = "\t")
data.table::fwrite(membership, file.path(out_dir, "cross_cohort_module_membership.tsv"), sep = "\t")
data.table::fwrite(membership_concordance, file.path(out_dir, "module_membership_concordance.tsv"), sep = "\t")
data.table::fwrite(clinical_assoc, file.path(out_dir, "adult_module_clinical_associations.tsv"), sep = "\t")
saveRDS(list(preservation = pres, test_eigengenes = test_mes, membership = membership), file.path(out_dir, "module_preservation_results.rds"))

module_focus <- trimws(strsplit(cfg$reference_modules, ",", fixed = TRUE)[[1]])
plot_pres <- preservation[module %in% module_focus]
p1 <- ggplot2::ggplot(plot_pres, ggplot2::aes(x = reorder(module, get(z_col)), y = get(z_col), fill = preservation_class)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_hline(yintercept = c(as.numeric(cfg$zsummary_moderate), as.numeric(cfg$zsummary_strong)), linetype = c(3, 2), color = "#555555") +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c(strong = "#1B9E77", moderate = "#E6AB02", weak_or_none = "#BDBDBD")) +
  ggplot2::labs(x = NULL, y = "Zsummary preservation", fill = NULL, title = "Reference-module preservation in adult SLE") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top")
ggplot2::ggsave(file.path(fig_dir, "fig02A_module_preservation.png"), p1, width = 7.4, height = 4.8, dpi = 320)

plot_mem <- membership[module %in% module_focus & is.finite(child_kME) & is.finite(adult_kME)]
p2 <- ggplot2::ggplot(plot_mem, ggplot2::aes(child_kME, adult_kME, color = module)) +
  ggplot2::geom_point(alpha = 0.22, size = 0.75) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#222222", linewidth = 0.6) +
  ggplot2::facet_wrap(~module, scales = "free") +
  ggplot2::labs(x = "Child cohort kME", y = "Adult cohort kME", title = "Cross-cohort module membership concordance") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "none")
ggplot2::ggsave(file.path(fig_dir, "fig02B_kME_concordance.png"), p2, width = 10.5, height = 4.2, dpi = 320)

if (nrow(clinical_assoc)) {
  p3 <- ggplot2::ggplot(clinical_assoc[module %in% module_focus], ggplot2::aes(trait, module, size = -log10(pmax(fdr, 1e-12)), color = beta)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
    ggplot2::labs(x = NULL, y = NULL, size = "-log10 FDR", color = "Beta", title = "Adult module-clinical associations") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  ggplot2::ggsave(file.path(fig_dir, "fig02C_module_clinical_bubble.png"), p3, width = 7.6, height = 4.6, dpi = 320)
}
