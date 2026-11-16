if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
project_root <- Sys.getenv("SLE_PROJECT_ROOT")
stopifnot(nzchar(project_root))
results_root <- file.path(project_root, "results", "3.2")
reference_root <- file.path(project_root, "reference", "3.2")
summary_root <- file.path(results_root, "04_summary")
dir.create(summary_root, recursive = TRUE, showWarnings = FALSE)

required <- c("data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

audit <- data.table::fread(file.path(results_root, "01_expression", "expression_audit.tsv"))
deg <- data.table::fread(file.path(results_root, "02_deg_gsea", "full_transcriptome_limma.tsv"))
gsea <- data.table::fread(file.path(results_root, "02_deg_gsea", "fgsea_all_collections.tsv"))
sensitivity <- data.table::fread(file.path(results_root, "02_deg_gsea", "covariate_sensitivity_summary.tsv"))
gsea_sensitivity <- data.table::fread(file.path(results_root, "02_deg_gsea", "fgsea_covariate_sensitivity.tsv"))
pres <- data.table::fread(file.path(results_root, "03_preservation", "module_preservation_statistics.tsv"))
membership <- data.table::fread(file.path(results_root, "03_preservation", "module_membership_concordance.tsv"))
clinical <- data.table::fread(file.path(results_root, "03_preservation", "adult_module_clinical_associations.tsv"))

z_col <- grep("^Zsummary", names(pres), value = TRUE)[1]
if (is.na(z_col)) stop("Current preservation result has no Zsummary column.")
get_frozen <- function(program_label) {
  gsea[gsea[["pathway"]] == paste0("FROZEN::", program_label)][1]
}
frozen_labels <- c("MAGENTA_IFN_ANTIVIRAL", "CYAN_NEUTROPHIL_IMMATURE", "GREENYELLOW_B_CELL_COMPLEMENT")
frozen <- data.table::rbindlist(lapply(frozen_labels, function(label) {
  x <- get_frozen(label)
  data.table::data.table(
    feature = label,
    NES = x$NES,
    FDR = x$padj,
    size = x$size
  )
}), fill = TRUE)

primary_count <- sum(deg$direction != "not_primary")
up_count <- sum(deg$direction == "SLE_up")
down_count <- sum(deg$direction == "SLE_down")
findings <- data.table::rbindlist(list(
  data.table::data.table(finding_id = "A01", layer = "input", feature = "GSE88884", estimate = sprintf("%d genes x %d samples", audit[metric == "selected_genes", value], audit[metric == "raw_samples", value]), uncertainty = "1,760 SLE; 60 healthy", interpretation = "adult full-transcriptome input"),
  data.table::data.table(finding_id = "A02", layer = "DEG", feature = "primary threshold genes", estimate = sprintf("%d total; %d up; %d down", primary_count, up_count, down_count), uncertainty = "BH FDR <= 0.05; abs(log2FC) >= 0.25", interpretation = "primary adult case-control signal"),
  frozen[, .(finding_id = paste0("A0", seq_len(.N) + 2L), layer = "frozen_program_GSEA", feature, estimate = sprintf("NES=%.6f", NES), uncertainty = sprintf("FDR=%.6g; n=%d", FDR, size), interpretation = ifelse(feature == "GREENYELLOW_B_CELL_COMPLEMENT", "not a positive adult case-control axis", "positive adult SLE enrichment"))],
  sensitivity[, .(finding_id = paste0("A1", seq_len(.N)), layer = "case_control_sensitivity", feature = model, estimate = sprintf("primary FDR replication=%.4f", primary_gene_fdr_replication), uncertainty = sprintf("n=%d SLE; n=%d healthy", n_sle, n_healthy), interpretation = "internal sensitivity; controls reused or female-restricted")],
  pres[module %in% c("magenta", "cyan", "greenyellow"), .(finding_id = paste0("P0", seq_len(.N)), layer = "module_preservation", feature = module, estimate = sprintf("Zsummary=%.6f", get(z_col)), uncertainty = sprintf("module size in common genes=%d; class=%s", moduleSize.z, preservation_class), interpretation = "adult SLE-only network preservation")],
  membership[module %in% c("magenta", "cyan", "greenyellow"), .(finding_id = paste0("P1", seq_len(.N)), layer = "kME_concordance", feature = module, estimate = sprintf("Spearman=%.6f", spearman_kME), uncertainty = sprintf("n genes=%d", n_genes), interpretation = "child-adult membership concordance")],
  clinical[module == "magenta" & trait == "C3", .(finding_id = "P20", layer = "adult_module_trait", feature = "magenta IFN vs C3", estimate = sprintf("beta=%.6f", beta), uncertainty = sprintf("n=%d; FDR=%.6g", n, fdr), interpretation = "unadjusted cross-sectional coherence signal")]
), fill = TRUE)
data.table::fwrite(findings, file.path(summary_root, "current_key_findings.tsv"), sep = "\t")

pre_path <- file.path(reference_root, "results", "module_preservation_pre_bicor", "module_preservation_statistics.tsv")
pre <- data.table::fread(pre_path)
pre_z <- grep("^Zsummary", names(pre), value = TRUE)[1]
comparison <- merge(
  pres[module %in% c("magenta", "cyan", "greenyellow"), .(module, current_zsummary = get(z_col), current_class = preservation_class)],
  pre[module %in% c("magenta", "cyan", "greenyellow"), .(module, pre_bicor_zsummary = get(pre_z), pre_bicor_class = preservation_class)],
  by = "module", all = TRUE
)
comparison[, delta_zsummary := current_zsummary - pre_bicor_zsummary]
data.table::fwrite(comparison, file.path(summary_root, "module_preservation_pre_post_comparison.tsv"), sep = "\t")

gsea_sensitivity_frozen <- gsea_sensitivity[collection == "FROZEN" & label %in% frozen_labels, .(model, label, NES, padj)]
data.table::fwrite(gsea_sensitivity_frozen, file.path(summary_root, "frozen_program_sensitivity.tsv"), sep = "\t")
cat("3.2 summary written.\n")
