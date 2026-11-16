if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
project_root <- Sys.getenv("SLE_PROJECT_ROOT")
results_root <- file.path(project_root, "results", "3.3")
grouping_root <- file.path(results_root, "01_grouping")
summary_root <- file.path(results_root, "02_summary")
dir.create(summary_root, recursive = TRUE, showWarnings = FALSE)

metrics <- data.table::fread(file.path(grouping_root, "cluster_number_stability.tsv"))
assignments <- data.table::fread(file.path(grouping_root, "baseline_subtype_assignments.tsv"))
assignment_audit <- data.table::fread(file.path(grouping_root, "validation_assignment_audit.tsv"))
method_audit <- data.table::fread(file.path(grouping_root, "distance_vs_correlation_assignment_audit.tsv"))
overlap <- data.table::fread(file.path(grouping_root, "trial_baseline_overlap_audit.tsv"))

selected <- metrics[order(selection_score, PAC, -mean_silhouette)][1]
group_counts <- assignments[, .N, by = .(cohort_partition, subtype)]
data.table::fwrite(group_counts, file.path(summary_root, "group_counts_by_partition.tsv"), sep = "\t")

findings <- data.table::rbindlist(list(
  data.table::data.table(
    result_id = c("3.3_K1", "3.3_K2", "3.3_K3"),
    analysis = "consensus_grouping",
    metric = c("selected_k", "PAC_at_selected_k", "mean_silhouette_at_selected_k"),
    value = c(selected$k, selected$PAC, selected$mean_silhouette),
    interpretation = c("working four-group solution", "lowest PAC", "incomplete separation")
  ),
  assignment_audit[cohort_partition != "discovery", data.table::data.table(
    result_id = ifelse(cohort_partition == "ILLUMINATE1_validation", "3.3_K4", "3.3_K5"),
    analysis = "heldout_assignment",
    metric = "low_distance_margin_fraction",
    value = low_confidence_fraction,
    interpretation = "held-out trial-linked assignment, not independent external replication"
  )],
  method_audit[, data.table::data.table(
    result_id = ifelse(cohort_partition == "ILLUMINATE1_validation", "3.3_K6", "3.3_K7"),
    analysis = "method_sensitivity",
    metric = "distance_vs_correlation_changed_fraction",
    value = changed_fraction,
    interpretation = "algorithm descriptions are not interchangeable"
  )]
), use.names = TRUE)
data.table::fwrite(findings, file.path(summary_root, "current_key_findings.tsv"), sep = "\t")

boundary <- data.table::data.table(
  topic = c("grouping ontology", "held-out samples", "assignment margin", "clinical comparisons", "downstream use"),
  supported_statement = c(
    "four partially separated program-based groupings",
    "trial-linked baseline samples can be assigned without refitting",
    "relative squared-distance margin is a heuristic assignment margin",
    "combined and partition-specific associations are descriptive",
    "frozen Euclidean assignments may be used by Results 3.4"
  ),
  unsupported_statement = c(
    "four independently reproduced biological disease states",
    "independent external cohort validation",
    "calibrated probability or confidence",
    "covariate-independent clinical effects",
    "validated treatment-response prediction"
  )
)
data.table::fwrite(boundary, file.path(summary_root, "interpretive_boundaries.tsv"), sep = "\t")

overlap_summary <- overlap[, .(
  matched_subjects = max(matched_subjects),
  minimum_cel_match_fraction = min(cel_match_fraction),
  minimum_program_score_pearson = min(program_score_pearson)
), by = trial]
data.table::fwrite(overlap_summary, file.path(summary_root, "trial_overlap_summary.tsv"), sep = "\t")
cat("3.3 summary tables written.\n")

