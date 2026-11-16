if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "3.1")
summary_dir <- file.path(results_root, "08_summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required for 3.1 result summaries.")
}

read_tsv <- function(path) data.table::fread(path, sep = "\t", na.strings = "NA")
write_tsv <- function(x, path) data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")

audit <- read_tsv(file.path(results_root, "07_module_selection_audit", "Table_S2_module_selection_audit.tsv"))
traits <- read_tsv(file.path(results_root, "06_wgcna", "module_trait_associations.tsv"))
stability <- read_tsv(file.path(results_root, "06_wgcna", "module_stability.tsv"))
wgcna_audit <- read_tsv(file.path(results_root, "06_wgcna", "wgcna_input_audit.tsv"))

selected_rows <- audit[audit[["selected"]] == TRUE][order(selection_rank)]
selected_modules <- selected_rows$module
pick <- function(module, family) {
  module_value <- module
  family_value <- family
  row <- traits[traits[["module"]] == module_value & traits[["model_family"]] == family_value][1]
  if (nrow(row) != 1L) stop("Missing module-trait row for ", module, " / ", family)
  row
}

module_rows <- lapply(selected_modules, function(module) {
  base <- pick(module, "activity_base")
  adjusted <- pick(module, "activity_composition_adjusted")
  c4 <- pick(module, "c4_composition_adjusted")
  stable <- stability[stability[["module"]] == module][1]
  selected_row <- selected_rows[selected_rows[["module"]] == module][1]
  data.table::data.table(
    module = module,
    biological_label = selected_row$biological_label,
    n_genes = selected_row$n_genes,
    activity_base_estimate = base$estimate,
    activity_base_fdr = base$fdr,
    activity_composition_adjusted_estimate = adjusted$estimate,
    activity_composition_adjusted_fdr = adjusted$fdr,
    c4_composition_adjusted_estimate = c4$estimate,
    c4_composition_adjusted_fdr = c4$fdr,
    median_kme_rank_correlation = stable$median_kme_rank_correlation,
    selected = TRUE
  )
})
write_tsv(data.table::rbindlist(module_rows), file.path(summary_dir, "selected_module_summary.tsv"))

metric <- function(name) wgcna_audit[metric == name, value][1]
findings <- data.table::data.table(
  finding = c(
    "analysis_population", "network", "selected_modules", "base_activity_association",
    "composition_adjusted_activity", "greenyellow_support", "stability", "reproduction"
  ),
  current_result = c(
    paste0(metric("input_visits"), " visits from ", metric("input_patients"), " pediatric SLE patients; ", metric("input_genes"), " input genes"),
    paste0(metric("selected_genes"), " high-variance genes; signed bicorrelation network; power ", metric("chosen_power"), "; ", metric("non_grey_modules"), " non-grey modules"),
    paste0(selected_rows[1, module], " ", selected_rows[1, n_genes], " genes; ", selected_rows[2, module], " ", selected_rows[2, n_genes], " genes; ", selected_rows[3, module], " ", selected_rows[3, n_genes], " genes"),
    "All three selected modules were positively associated with within-patient SLEDAI at FDR < 0.05",
    paste0("magenta and cyan remained significant; greenyellow had adjusted SLEDAI FDR=", selected_rows[selected_rows[["module"]] == "greenyellow", minimum_adjusted_fdr][1]),
    "greenyellow retained an inverse within-patient C4 association after composition adjustment, FDR=8.79e-07",
    paste0("Median kME rank correlations were ", paste(format(stability[module %in% selected_modules, median_kme_rank_correlation], digits = 4), collapse = ", "), "; all 25 bootstrap correlations were at least 0.70"),
    "Current WGCNA core tables numerically match the reference snapshot"
  ),
  evidence_file = c(
    "../06_wgcna/wgcna_input_audit.tsv", "../06_wgcna/wgcna_input_audit.tsv",
    "../07_module_selection_audit/Table_S2_module_selection_audit.tsv", "selected_module_summary.tsv",
    "selected_module_summary.tsv", "selected_module_summary.tsv", "../06_wgcna/module_stability.tsv",
    "../../../analysis/verify_3_1_reproduction.R"
  )
)
write_tsv(findings, file.path(summary_dir, "current_key_findings.tsv"))

write_tsv(data.table::data.table(
  topic = c("module_selection", "bootstrap", "greenyellow", "low_c4", "downstream_impact"),
  boundary = c(
    "The three-module choice was an exploratory discovery-stage selection, not a pre-WGCNA preregistration.",
    "The 25 patient-level bootstraps recalculated kME within fixed modules; they did not rebuild and rematch the full network.",
    "Do not state that greenyellow remained significantly associated with SLEDAI after composition adjustment; its adjusted SLEDAI FDR was 0.0526.",
    "The adjusted inverse C4 association supports greenyellow selection but does not establish causality or a cell-intrinsic mechanism.",
    "Module membership and sizes were unchanged, so downstream projections are numerically unchanged; only the strength and wording of the evidence claim changes."
  )
), file.path(summary_dir, "interpretive_boundaries.tsv"))

write_tsv(data.table::data.table(
  category = c("input_audit", "longitudinal_models", "wgcna_input", "module_associations", "module_stability", "module_membership", "selection_audit", "revised_text", "interpretation"),
  path = c(
    "../00_input_audit/input_audit_summary.tsv", "../03_longitudinal/longitudinal_model_results.tsv",
    "../06_wgcna/wgcna_input_audit.tsv", "../06_wgcna/module_trait_associations.tsv",
    "../06_wgcna/module_stability.tsv", "../06_wgcna/module_membership.tsv",
    "../07_module_selection_audit/Table_S2_module_selection_audit.tsv",
    "../../../manuscript/revised_sections/03_Results_3.1_revised.md",
    "../../../docs/3.1_modification_and_interpretation.md"
  ),
  description = c(
    "Included-array and patient audit", "Program-level longitudinal models", "WGCNA sample, gene, power, and bootstrap counts",
    "Base and composition-adjusted module-trait models", "Patient-bootstrap kME rank-correlation stability",
    "Gene-module membership and kME", "Selection gates, ranking, labels, and reasons",
    "Revised manuscript Results 3.1 text", "Chinese modification and interpretation record"
  )
), file.path(summary_dir, "result_manifest.tsv"))

message("3.1 summary outputs written to ", summary_dir)
