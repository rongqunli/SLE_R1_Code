if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
project_root <- Sys.getenv("SLE_PROJECT_ROOT")
analysis_root <- file.path(project_root, "pipeline", "3.3")
results_root <- file.path(project_root, "results", "3.3")
cfg_tab <- data.table::fread(file.path(analysis_root, "config", "analysis_contract.tsv"))
cfg <- setNames(as.list(cfg_tab$value), cfg_tab$key)
resolve <- function(path) normalizePath(file.path(analysis_root, path), winslash = "/", mustWork = TRUE)
set.seed(as.integer(cfg$random_seed))

base <- readRDS(resolve(cfg$gse88884_object))
trial1 <- readRDS(resolve(cfg$gse88885_object))
trial2 <- readRDS(resolve(cfg$gse88886_object))
clinical <- data.table::fread(
  resolve(cfg$frozen_metadata), sep = "\t",
  na.strings = c("", "NA", "not_available", "not_applicable", "--")
)

scores <- data.table::as.data.table(base$program_scores)
roles <- trimws(strsplit(cfg$feature_roles, ",", fixed = TRUE)[[1]])
scores <- scores[condition == "SLE" & time == "baseline" & analysis_role %in% roles]
wide <- data.table::dcast(scores, sample_id + subject_id ~ program, value.var = "score_z")
feature_cols <- setdiff(names(wide), c("sample_id", "subject_id"))
if (length(feature_cols) != 3L) stop("Expected exactly three primary frozen programs")
wide <- wide[stats::complete.cases(wide[, ..feature_cols])]

trial1_ids <- unique(data.table::as.data.table(trial1$metadata)$subject_id)
trial2_ids <- unique(data.table::as.data.table(trial2$metadata)$subject_id)
wide[, cohort_partition := data.table::fcase(
  subject_id %in% trial1_ids, "ILLUMINATE1_validation",
  subject_id %in% trial2_ids, "ILLUMINATE2_validation",
  default = "discovery"
)]
if (anyDuplicated(wide$subject_id)) stop("Baseline grouping table contains duplicated subjects")

discovery <- wide[cohort_partition == "discovery"]
x_raw <- as.matrix(discovery[, ..feature_cols])
center <- colMeans(x_raw)
scale_value <- apply(x_raw, 2, stats::sd)
if (any(!is.finite(scale_value) | scale_value == 0)) stop("Invalid grouping feature scale")
x <- sweep(sweep(x_raw, 2, center, "-"), 2, scale_value, "/")
rownames(x) <- discovery$subject_id

consensus_kmeans <- function(x, k, n_resamples, sample_fraction, seed) {
  set.seed(seed + k)
  n <- nrow(x)
  together <- matrix(0L, n, n)
  cocluster <- matrix(0L, n, n)
  for (b in seq_len(n_resamples)) {
    ids <- sort(sample.int(n, max(k * 5L, floor(n * sample_fraction)), replace = FALSE))
    km <- stats::kmeans(x[ids, , drop = FALSE], centers = k, nstart = 30, iter.max = 100)
    together[ids, ids] <- together[ids, ids] + 1L
    for (cl in seq_len(k)) {
      in_cl <- ids[km$cluster == cl]
      cocluster[in_cl, in_cl] <- cocluster[in_cl, in_cl] + 1L
    }
  }
  consensus <- cocluster / pmax(together, 1L)
  dimnames(consensus) <- list(rownames(x), rownames(x))
  diag(consensus) <- 1
  hc <- stats::hclust(stats::as.dist(1 - consensus), method = "average")
  cluster_id <- stats::cutree(hc, k = k)
  sil <- cluster::silhouette(cluster_id, stats::dist(x))
  off_diag <- consensus[lower.tri(consensus)]
  list(
    k = k,
    consensus = consensus,
    hclust = hc,
    cluster = cluster_id,
    pac = mean(off_diag > as.numeric(cfg$pac_lower) & off_diag < as.numeric(cfg$pac_upper), na.rm = TRUE),
    silhouette = mean(sil[, "sil_width"], na.rm = TRUE),
    cdf_area = mean(off_diag, na.rm = TRUE)
  )
}

k_values <- seq.int(as.integer(cfg$k_min), as.integer(cfg$k_max))
solutions <- lapply(k_values, function(k) consensus_kmeans(
  x, k, as.integer(cfg$consensus_resamples), as.numeric(cfg$sample_fraction), as.integer(cfg$random_seed)
))
metrics <- data.table::rbindlist(lapply(solutions, function(z) data.table::data.table(
  k = z$k, PAC = z$pac, mean_silhouette = z$silhouette, consensus_area = z$cdf_area
)))
metrics[, selection_score := data.table::frank(PAC, ties.method = "average") + data.table::frank(-mean_silhouette, ties.method = "average")]
eligible <- metrics[mean_silhouette >= as.numeric(cfg$minimum_silhouette)]
if (!nrow(eligible)) warning("No k reached the silhouette threshold")
selected_k <- (if (nrow(eligible)) eligible else metrics)[order(selection_score, PAC, -mean_silhouette), k][1]
solution <- solutions[[match(selected_k, k_values)]]

discovery[, raw_cluster := solution$cluster[match(subject_id, names(solution$cluster))]]
centroids <- data.table::as.data.table(x)[, raw_cluster := discovery$raw_cluster][, lapply(.SD, mean), by = raw_cluster, .SDcols = feature_cols]
dominant <- apply(as.matrix(centroids[, ..feature_cols]), 1, function(v) feature_cols[which.max(v)])
order_idx <- order(dominant, centroids$raw_cluster)
centroids <- centroids[order_idx]
dominant <- dominant[order_idx]
program_label <- c(
  CYAN_NEUTROPHIL_IMMATURE = "NEUTROPHIL_HIGH",
  GREENYELLOW_B_CELL_COMPLEMENT = "B_CELL_HIGH",
  MAGENTA_IFN_ANTIVIRAL = "IFN_HIGH"
)
centroid_matrix <- as.matrix(centroids[, ..feature_cols])
biological_label <- ifelse(
  apply(centroid_matrix, 1, max) < 0,
  "IMMUNE_LOW",
  unname(program_label[dominant])
)
centroids[, subtype := sprintf("S%02d_%s", seq_len(.N), biological_label)]
cluster_label <- setNames(centroids$subtype, centroids$raw_cluster)
discovery[, subtype := unname(cluster_label[as.character(raw_cluster)])]

nearest_centroid <- function(tab, centroids, center, scale_value) {
  m <- as.matrix(tab[, ..feature_cols])
  m <- sweep(sweep(m, 2, center, "-"), 2, scale_value, "/")
  cmat <- as.matrix(centroids[, ..feature_cols])
  d <- sapply(seq_len(nrow(cmat)), function(j) rowSums((m - matrix(cmat[j, ], nrow(m), ncol(m), byrow = TRUE))^2))
  if (is.null(dim(d))) d <- matrix(d, ncol = nrow(cmat))
  best <- max.col(-d, ties.method = "first")
  sorted <- t(apply(d, 1, sort))
  margin <- if (ncol(sorted) >= 2L) (sorted[, 2] - sorted[, 1]) / pmax(sorted[, 2], 1e-8) else 1
  data.table::data.table(
    subject_id = tab$subject_id,
    subtype = centroids$subtype[best],
    centroid_distance = d[cbind(seq_len(nrow(d)), best)],
    assignment_margin = margin
  )
}

validation <- wide[cohort_partition != "discovery"]
mapped <- nearest_centroid(validation, centroids, center, scale_value)
validation <- merge(validation, mapped, by = "subject_id", all.x = TRUE)
discovery_margin <- nearest_centroid(discovery, centroids, center, scale_value)
discovery[, `:=`(centroid_distance = 0, assignment_margin = 1)]
assignments <- data.table::rbindlist(list(
  discovery[, c("sample_id", "subject_id", "cohort_partition", feature_cols, "subtype", "centroid_distance", "assignment_margin"), with = FALSE],
  validation[, c("sample_id", "subject_id", "cohort_partition", feature_cols, "subtype", "centroid_distance", "assignment_margin"), with = FALSE]
), use.names = TRUE)
assignments[, low_confidence := assignment_margin < as.numeric(cfg$minimum_assignment_margin)]

clinical_keep <- clinical[, .(sample_id, patient_id, SLEDAI, C3, C4, anti_dsDNA, age, sex, race, batch, cohort)]
annotated <- merge(assignments, clinical_keep, by = "sample_id", all.x = TRUE)
numeric_traits <- c("SLEDAI", "C3", "C4", "anti_dsDNA", "age")
clinical_tests <- data.table::rbindlist(lapply(numeric_traits, function(trait) {
  dat <- annotated[, .(subtype, value = suppressWarnings(as.numeric(get(trait))))]
  dat <- dat[is.finite(value)]
  kt <- stats::kruskal.test(value ~ subtype, data = dat)
  data.table::data.table(trait = trait, test = "Kruskal-Wallis", n = nrow(dat), statistic = unname(kt$statistic), p_value = kt$p.value)
}), fill = TRUE)
categorical_traits <- c("sex", "race", "batch")
clinical_tests <- data.table::rbindlist(c(list(clinical_tests), lapply(categorical_traits, function(trait) {
  dat <- annotated[!is.na(get(trait)) & nzchar(get(trait)), .(subtype, value = get(trait))]
  tab <- table(dat$subtype, dat$value)
  ct <- suppressWarnings(stats::chisq.test(tab, simulate.p.value = TRUE, B = 5000))
  cramer_v <- sqrt(unname(ct$statistic) / (sum(tab) * min(nrow(tab) - 1L, ncol(tab) - 1L)))
  data.table::data.table(
    trait = trait, test = "Chi-square simulated", n = sum(tab),
    statistic = unname(ct$statistic), p_value = ct$p.value,
    effect_size = cramer_v, effect_size_metric = "Cramer's V"
  )
})), fill = TRUE)
clinical_tests[, fdr := stats::p.adjust(p_value, method = "BH")]
assignment_audit <- annotated[, .(
  subjects = .N,
  low_confidence_subjects = sum(low_confidence),
  low_confidence_fraction = mean(low_confidence),
  median_assignment_margin = stats::median(assignment_margin),
  q10_assignment_margin = stats::quantile(assignment_margin, 0.10)
), by = cohort_partition]

out_dir <- file.path(results_root, "01_grouping")
fig_dir <- file.path(results_root, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(metrics, file.path(out_dir, "cluster_number_stability.tsv"), sep = "\t")
data.table::fwrite(centroids, file.path(out_dir, "frozen_subtype_centroids.tsv"), sep = "\t")
data.table::fwrite(annotated, file.path(out_dir, "baseline_subtype_assignments.tsv"), sep = "\t")
data.table::fwrite(clinical_tests, file.path(out_dir, "post_freeze_clinical_associations.tsv"), sep = "\t")
data.table::fwrite(assignment_audit, file.path(out_dir, "validation_assignment_audit.tsv"), sep = "\t")

group_assignments <- data.table::copy(annotated)
data.table::setnames(group_assignments, c("subtype", "low_confidence"), c("grouping", "low_margin"))
data.table::fwrite(group_assignments, file.path(out_dir, "baseline_group_assignments.tsv"), sep = "\t")
group_centroids <- data.table::copy(centroids)
data.table::setnames(group_centroids, "subtype", "grouping")
data.table::fwrite(group_centroids, file.path(out_dir, "frozen_group_centroids.tsv"), sep = "\t")

discovery_margin_audit <- merge(
  discovery[, .(subject_id, consensus_grouping = subtype)],
  discovery_margin[, .(subject_id, centroid_grouping = subtype, centroid_distance, assignment_margin)],
  by = "subject_id"
)
discovery_margin_audit[, grouping_agrees := consensus_grouping == centroid_grouping]
data.table::fwrite(discovery_margin_audit, file.path(out_dir, "discovery_centroid_margin_audit.tsv"), sep = "\t")

scaled_all <- sweep(sweep(as.matrix(wide[, ..feature_cols]), 2, center, "-"), 2, scale_value, "/")
rownames(scaled_all) <- wide$subject_id
fixed_silhouette <- data.table::rbindlist(lapply(unique(annotated$cohort_partition), function(partition) {
  ids <- annotated[cohort_partition == partition, subject_id]
  labels <- annotated[cohort_partition == partition, subtype]
  sil <- cluster::silhouette(as.integer(factor(labels)), stats::dist(scaled_all[ids, , drop = FALSE]))
  data.table::data.table(
    cohort_partition = partition,
    n = length(ids),
    mean_silhouette = mean(sil[, "sil_width"]),
    negative_silhouette_fraction = mean(sil[, "sil_width"] < 0)
  )
}))
data.table::fwrite(fixed_silhouette, file.path(out_dir, "fixed_grouping_silhouette_by_partition.tsv"), sep = "\t")

de_novo <- data.table::rbindlist(lapply(c("ILLUMINATE1_validation", "ILLUMINATE2_validation"), function(partition) {
  ids <- annotated[cohort_partition == partition, subject_id]
  mat <- scaled_all[ids, , drop = FALSE]
  data.table::rbindlist(lapply(k_values, function(k) {
    set.seed(as.integer(cfg$random_seed) + k)
    km <- stats::kmeans(mat, centers = k, nstart = 100, iter.max = 100)
    sil <- cluster::silhouette(km$cluster, stats::dist(mat))
    data.table::data.table(cohort_partition = partition, k = k, mean_silhouette = mean(sil[, "sil_width"]))
  }))
}))
data.table::fwrite(de_novo, file.path(out_dir, "heldout_de_novo_kmeans_silhouette.tsv"), sep = "\t")

clinical_by_partition <- data.table::rbindlist(lapply(unique(annotated$cohort_partition), function(partition) {
  result <- data.table::rbindlist(lapply(numeric_traits, function(trait) {
    dat <- annotated[cohort_partition == partition, .(subtype, value = suppressWarnings(as.numeric(get(trait))))]
    dat <- dat[is.finite(value)]
    kt <- stats::kruskal.test(value ~ subtype, data = dat)
    groups <- data.table::uniqueN(dat$subtype)
    data.table::data.table(
      cohort_partition = partition, trait = trait, n = nrow(dat),
      statistic = unname(kt$statistic), p_value = kt$p.value,
      epsilon_squared = (unname(kt$statistic) - groups + 1) / (nrow(dat) - groups)
    )
  }))
  result[, fdr_within_partition := stats::p.adjust(p_value, method = "BH")]
  result
}))
data.table::fwrite(clinical_by_partition, file.path(out_dir, "clinical_associations_by_partition.tsv"), sep = "\t")

trial_overlap <- function(trial_object, trial_name) {
  trial_scores <- data.table::as.data.table(trial_object$program_scores)[
    condition == "SLE" & time == "baseline" & analysis_role %in% roles,
    .(subject_id, program, trial_sample_id = sample_id, trial_expression_column = expression_column, trial_score = score_z)
  ]
  base_scores <- scores[, .(subject_id, program, base_sample_id = sample_id, base_expression_column = expression_column, base_score = score_z)]
  joined <- merge(base_scores, trial_scores, by = c("subject_id", "program"))
  joined[, .(
    trial = trial_name,
    matched_subjects = data.table::uniqueN(subject_id),
    matched_score_rows = .N,
    exact_cel_matches = sum(base_expression_column == trial_expression_column),
    cel_match_fraction = mean(base_expression_column == trial_expression_column),
    same_geo_sample_id_fraction = mean(base_sample_id == trial_sample_id),
    program_score_pearson = stats::cor(base_score, trial_score)
  ), by = program]
}
overlap_audit <- data.table::rbindlist(list(
  trial_overlap(trial1, "ILLUMINATE-1"),
  trial_overlap(trial2, "ILLUMINATE-2")
))
data.table::fwrite(overlap_audit, file.path(out_dir, "trial_baseline_overlap_audit.tsv"), sep = "\t")

validation_scaled <- scaled_all[validation$subject_id, , drop = FALSE]
cor_matrix <- t(apply(validation_scaled, 1, function(v) apply(centroid_matrix, 1, function(z) stats::cor(v, z))))
cor_best <- max.col(cor_matrix, ties.method = "first")
cor_sorted <- t(apply(cor_matrix, 1, sort, decreasing = TRUE))
cor_assignments <- data.table::data.table(
  subject_id = validation$subject_id,
  cohort_partition = validation$cohort_partition,
  distance_grouping = validation$subtype,
  correlation_grouping = centroids$subtype[cor_best],
  correlation_margin = cor_sorted[, 1] - cor_sorted[, 2]
)
cor_assignments[, grouping_changed := distance_grouping != correlation_grouping]
data.table::fwrite(cor_assignments, file.path(out_dir, "alternative_correlation_assignments.tsv"), sep = "\t")
method_audit <- cor_assignments[, .(
  n = .N,
  changed_assignments = sum(grouping_changed),
  changed_fraction = mean(grouping_changed),
  correlation_low_margin_n = sum(correlation_margin < as.numeric(cfg$minimum_assignment_margin)),
  correlation_low_margin_fraction = mean(correlation_margin < as.numeric(cfg$minimum_assignment_margin)),
  median_correlation_margin = stats::median(correlation_margin)
), by = cohort_partition]
method_audit <- merge(
  method_audit,
  assignment_audit[cohort_partition != "discovery", .(
    cohort_partition,
    distance_low_margin_n = low_confidence_subjects,
    distance_low_margin_fraction = low_confidence_fraction,
    median_distance_margin = median_assignment_margin
  )],
  by = "cohort_partition"
)
data.table::fwrite(method_audit, file.path(out_dir, "distance_vs_correlation_assignment_audit.tsv"), sep = "\t")

method_definition <- data.table::data.table(
  item = c("distance", "assignment", "margin", "low_margin_threshold", "interpretation"),
  implementation = c(
    "squared Euclidean distance in discovery-standardized three-program space",
    "nearest frozen centroid",
    "(d2-d1)/d2, where d1 and d2 are the smallest and second-smallest squared distances",
    cfg$minimum_assignment_margin,
    "held-out trial-linked assignment; not independent external replication"
  )
)
data.table::fwrite(method_definition, file.path(out_dir, "assignment_method_definition.tsv"), sep = "\t")

saveRDS(list(
  selected_k = selected_k, metrics = metrics, consensus = solution$consensus,
  hclust = solution$hclust, centroids = centroids, center = center, scale = scale_value,
  features = feature_cols, assignments = annotated
), file.path(out_dir, "frozen_grouping_model.rds"))

metric_long <- data.table::melt(metrics, id.vars = "k", measure.vars = c("PAC", "mean_silhouette"), variable.name = "metric")
p1 <- ggplot2::ggplot(metric_long, ggplot2::aes(k, value, color = metric)) +
  ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2.4) +
  ggplot2::geom_vline(xintercept = selected_k, linetype = 2, color = "#333333") +
  ggplot2::scale_color_manual(values = c(PAC = "#D95F02", mean_silhouette = "#1B9E77")) +
  ggplot2::scale_x_continuous(breaks = k_values) +
  ggplot2::labs(x = "Number of groupings", y = "Metric value", color = NULL, title = "Grouping stability without clinical-outcome selection") +
  ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(fig_dir, "fig03A_grouping_stability.png"), p1, width = 7.4, height = 4.8, dpi = 320)

ord <- solution$hclust$order
heat_mat <- t(x[ord, , drop = FALSE])
ha <- ComplexHeatmap::HeatmapAnnotation(
  grouping = discovery$subtype[ord],
  col = list(grouping = setNames(circlize::rand_color(length(unique(discovery$subtype)), luminosity = "bright"), sort(unique(discovery$subtype)))),
  show_annotation_name = FALSE
)
grDevices::png(file.path(fig_dir, "fig03B_grouping_heatmap.png"), width = 3600, height = 1700, res = 300)
ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
  heat_mat, name = "z", top_annotation = ha, cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE, col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
  column_title = paste0("Discovery SLE subjects; selected k = ", selected_k),
  row_names_gp = grid::gpar(fontsize = 9)
))
grDevices::dev.off()

pca <- stats::prcomp(x, center = FALSE, scale. = FALSE)
pca_df <- data.table::data.table(subject_id = rownames(x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], grouping = discovery$subtype)
p2 <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = grouping)) +
  ggplot2::geom_point(alpha = 0.6, size = 1.2) +
  ggplot2::stat_ellipse(linewidth = 0.7, level = 0.80) +
  ggplot2::labs(color = NULL, title = "Program-based groupings in three-program space") +
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
ggplot2::ggsave(file.path(fig_dir, "fig03C_grouping_embedding.png"), p2, width = 8.4, height = 6.1, dpi = 320)

alluvial <- annotated[, .N, by = .(cohort_partition, subtype, low_confidence)]
alluvial[, partition_label := data.table::fcase(
  cohort_partition == "discovery", "Discovery",
  cohort_partition == "ILLUMINATE1_validation", "ILLUMINATE-1 held out",
  cohort_partition == "ILLUMINATE2_validation", "ILLUMINATE-2 held out",
  default = cohort_partition
)]
alluvial[, grouping_label := sub("^S[0-9]+_", "", subtype)]
alluvial[, margin_label := ifelse(low_confidence, "Low margin", "Higher margin")]
p3 <- ggplot2::ggplot(alluvial, ggplot2::aes(axis1 = partition_label, axis2 = grouping_label, axis3 = margin_label, y = N)) +
  ggalluvial::geom_alluvium(ggplot2::aes(fill = subtype), width = 0.16, alpha = 0.76) +
  ggalluvial::geom_stratum(width = 0.18, fill = "#F2F2F2", color = "#666666") +
  ggalluvial::stat_stratum(geom = "text", ggplot2::aes(label = ggplot2::after_stat(stratum)), size = 3.1) +
  ggplot2::scale_x_discrete(limits = c("Partition", "Program grouping", "Assignment margin"), expand = c(.16, .16)) +
  ggplot2::labs(y = "Subjects", fill = NULL, title = "Discovery and held-out trial-linked assignments") +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(18, 70, 28, 70))
ggplot2::ggsave(file.path(fig_dir, "fig03D_heldout_assignment_alluvial.png"), p3, width = 12.8, height = 6.5, dpi = 320, bg = "white")

cat("3.3 grouping analysis completed with selected k=", selected_k, ".\n", sep = "")

