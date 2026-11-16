if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "3.1")
input_dir <- file.path(results_root, "06_wgcna")
output_dir <- file.path(results_root, "07_module_selection_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required for the module-selection audit.")
}

traits_file <- file.path(input_dir, "module_trait_associations.tsv")
stability_file <- file.path(input_dir, "module_stability.tsv")
membership_file <- file.path(input_dir, "module_membership.tsv")
required <- c(traits_file, stability_file, membership_file)
if (!all(file.exists(required))) {
  stop("Run pipeline step 06_wgcna first. Missing: ", paste(required[!file.exists(required)], collapse = ", "))
}

traits <- data.table::fread(traits_file, sep = "\t", na.strings = "NA")
traits <- traits[model_family %in% c(
  "activity_composition_adjusted", "renal_composition_adjusted",
  "c3_composition_adjusted", "c4_composition_adjusted"
)]
traits[, aligned_t := ifelse(grepl("^c[34]_", model_family), -t_statistic, t_statistic)]
traits[, axis := c(
  activity_composition_adjusted = "SLEDAI adjusted",
  renal_composition_adjusted = "Renal activity adjusted",
  c3_composition_adjusted = "Low C3 adjusted",
  c4_composition_adjusted = "Low C4 adjusted"
)[model_family]]

stability <- data.table::fread(stability_file, sep = "\t", na.strings = "NA")
membership <- data.table::fread(membership_file, sep = "\t", na.strings = "NA")
module_sizes <- membership[, .(module_n_genes = .N), by = module]

summary <- traits[, .(
  selection_score = max(aligned_t),
  best_axis = axis[which.max(aligned_t)],
  minimum_adjusted_fdr = min(fdr),
  significant_adjusted_axes = sum(fdr < 0.05)
), by = module]
summary <- merge(summary, stability, by = "module", all.x = TRUE)
summary <- merge(summary, module_sizes, by = "module", all.x = TRUE)
summary[, n_genes := module_n_genes]
summary[, eligible := significant_adjusted_axes >= 1 &
  median_kme_rank_correlation >= 0.97 & proportion_correlation_ge_0_7 == 1]
summary[, positive_candidate := eligible & selection_score > 0]
data.table::setorder(summary, -positive_candidate, -selection_score, module)
summary[, selection_rank := NA_integer_]
summary[positive_candidate == TRUE, selection_rank := seq_len(.N)]
summary[, selected := positive_candidate & selection_rank <= 3]
labels <- c(magenta = "IFN/antiviral", greenyellow = "B-cell/complement", cyan = "Immature-neutrophil")
summary[, biological_label := unname(labels[module])]
summary[, selection_reason := ifelse(
  selected,
  data.table::fifelse(module == "magenta",
          "Highest positive adjusted activity-aligned statistic; IFN/antiviral hubs",
          data.table::fifelse(module == "greenyellow",
                  "Highest positive adjusted low-C4 statistic; B-cell hubs and low-complement association",
                  "Third-ranked positive adjusted statistic; immature-neutrophil hubs and renal/activity support")),
  data.table::fifelse(positive_candidate,
          "Passed evidence gates but ranked below the top-three transport set",
          data.table::fifelse(significant_adjusted_axes >= 1,
                  "Passed evidence gates but lacked a positive activity-aligned association",
                  "Did not pass the adjusted-association and stability evidence gates"))
)]

out <- summary[, .(
  module, n_genes, selection_score, best_axis, minimum_adjusted_fdr,
  significant_adjusted_axes, median_kme_rank_correlation,
  q025_kme_rank_correlation, q975_kme_rank_correlation,
  proportion_correlation_ge_0_7, eligible, positive_candidate,
  selection_rank, selected, biological_label, selection_reason
)]
data.table::fwrite(out, file.path(output_dir, "Table_S2_module_selection_audit.tsv"), sep = "\t", quote = FALSE, na = "NA")

selected <- out[selected == TRUE]
if (!identical(selected$module, c("magenta", "greenyellow", "cyan"))) {
  stop("The audited selected modules differ from the manuscript v3 set.")
}
expected_sizes <- c(magenta = 157L, greenyellow = 64L, cyan = 35L)
if (!all(selected$n_genes[match(names(expected_sizes), selected$module)] == expected_sizes)) {
  stop("The audited module sizes differ from the manuscript v3 set.")
}

writeLines(c(
  "3.1 module-selection audit",
  paste0("Input WGCNA directory: ", input_dir),
  "Selection score: maximum activity-aligned t across four composition-adjusted families.",
  "Eligibility: at least one adjusted FDR < 0.05, median kME rank correlation >= 0.97, and all 25 bootstrap correlations >= 0.70.",
  "Selected modules: magenta, greenyellow, cyan.",
  "This script audits the exploratory discovery-stage choice; it does not alter module membership."
), file.path(output_dir, "selection_audit_summary.txt"))
message("Wrote ", file.path(output_dir, "Table_S2_module_selection_audit.tsv"))
