if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[[1]]), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(
  c("Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "patchwork"),
  "GSE135779 QC and clustering"
)

input_file <- file.path(cfg$paths$processed, "00_raw_seurat.rds")
output_file <- file.path(cfg$paths$processed, "01_qc_clustered_seurat.rds")
output_dir <- file.path(cfg$paths$results, "03_qc_cluster")
checkpoint_dir <- file.path(cfg$paths$processed, "checkpoints", "03_qc_cluster")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input object not found: ", input_file, ". Run 02_build_seurat_object.R first.")
}

signature_code_files <- c(
  file.path(PIPELINE_ROOT, "scripts", "03_qc_cluster.R"),
  file.path(PIPELINE_ROOT, "scripts", "03_recluster_from_pca.R"),
  file.path(PIPELINE_ROOT, "scripts", "_common.R"),
  file.path(PIPELINE_ROOT, "config", "config.R")
)
step_signature <- build_step_signature(
  input_files = input_file,
  extra = list(
    step = "03_qc_cluster",
    qc = cfg$qc,
    dimension_reduction = cfg$dimension_reduction
  ),
  code_files = signature_code_files
)
if (step_is_complete(output_dir, step_signature, output_file)) {
  log_message("Step 03 is already complete and its signature matches. Use --force to rebuild.")
  quit(save = "no", status = 0L)
}
initialize_step(
  step_signature,
  checkpoint_dir,
  owned_paths = c(output_dir, output_file)
)

set.seed(cfg$seed)

adaptive_upper <- function(x, hard_max = Inf, multiplier = 4) {
  x <- as.numeric(x)
  center <- stats::median(x, na.rm = TRUE)
  spread <- stats::mad(
    x,
    center = center,
    constant = cfg$qc$mad_constant,
    na.rm = TRUE
  )
  if (!is.finite(spread) || spread <= 0) {
    spread <- stats::IQR(x, na.rm = TRUE) / cfg$qc$iqr_normalizer
  }
  if (!is.finite(spread) || spread <= 0) {
    candidate <- stats::quantile(
      x,
      cfg$qc$fallback_upper_quantile,
      na.rm = TRUE,
      names = FALSE
    )
  } else {
    candidate <- center + multiplier * spread
  }
  min(as.numeric(candidate), hard_max)
}

basic_checkpoint <- file.path(checkpoint_dir, "01_after_basic_qc.rds")
basic_state <- load_checkpoint_rds(basic_checkpoint, step_signature)
if (is.null(basic_state)) {
  object <- readRDS(input_file)
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
  object$cohort_condition <- paste(object$cohort, object$condition, sep = "_")

  qc_before <- object@meta.data
  qc_before$cell_id <- rownames(qc_before)
  sample_ids <- sort(unique(as.character(object$sample_id)))
  thresholds <- lapply(sample_ids, function(sample_id) {
    md <- qc_before[qc_before$sample_id == sample_id, , drop = FALSE]
    data.frame(
      sample_id = sample_id,
      cells_before = nrow(md),
      nFeature_min = cfg$qc$min_features,
      nFeature_max = adaptive_upper(
        md$nFeature_RNA,
        hard_max = cfg$qc$max_features_hard,
        multiplier = cfg$qc$mad_multiplier
      ),
      nCount_max = adaptive_upper(
        md$nCount_RNA,
        hard_max = Inf,
        multiplier = cfg$qc$mad_multiplier
      ),
      percent_mt_max = cfg$qc$max_percent_mt,
      stringsAsFactors = FALSE
    )
  })
  thresholds <- do.call(rbind, thresholds)

  threshold_index <- match(object$sample_id, thresholds$sample_id)
  fail_low_feature <- object$nFeature_RNA < thresholds$nFeature_min[threshold_index]
  fail_high_feature <- object$nFeature_RNA > thresholds$nFeature_max[threshold_index]
  fail_high_count <- object$nCount_RNA > thresholds$nCount_max[threshold_index]
  fail_high_mt <- object$percent.mt > thresholds$percent_mt_max[threshold_index]
  failure_flags <- data.frame(
    sample_id = as.character(object$sample_id),
    low_feature = fail_low_feature,
    high_feature = fail_high_feature,
    high_count = fail_high_count,
    high_mt = fail_high_mt,
    stringsAsFactors = FALSE
  )
  flag_columns <- c("low_feature", "high_feature", "high_count", "high_mt")
  for (column in flag_columns) {
    failure_flags[[column]][is.na(failure_flags[[column]])] <- FALSE
  }
  failure_summary <- data.table::as.data.table(failure_flags)[, .(
    cells = .N,
    failed_low_feature = sum(low_feature),
    failed_high_feature = sum(high_feature),
    failed_high_count = sum(high_count),
    failed_high_mt = sum(high_mt),
    failed_any = sum(low_feature | high_feature | high_count | high_mt)
  ), by = sample_id]
  keep_qc <- !(fail_low_feature | fail_high_feature | fail_high_count | fail_high_mt)
  keep_qc[is.na(keep_qc)] <- FALSE
  object$basic_qc_pass <- keep_qc
  object <- subset(object, cells = colnames(object)[keep_qc])

  cells_after_basic <- table(as.character(object$sample_id))
  thresholds$cells_after_basic_qc <- as.integer(cells_after_basic[thresholds$sample_id])
  thresholds$cells_after_basic_qc[is.na(thresholds$cells_after_basic_qc)] <- 0L
  thresholds$basic_qc_retained_fraction <- thresholds$cells_after_basic_qc / thresholds$cells_before
  too_small <- thresholds$sample_id[
    thresholds$cells_after_basic_qc < cfg$qc$min_cells_per_sample_after_qc
  ]
  if (length(too_small) > 0L) {
    warning("Samples with few cells after basic QC: ", paste(too_small, collapse = ", "))
    if (isTRUE(cfg$qc$drop_samples_below_min_cells)) {
      keep_cells <- colnames(object)[!as.character(object$sample_id) %in% too_small]
      object <- subset(object, cells = keep_cells)
      thresholds$sample_retained <- !thresholds$sample_id %in% too_small
    }
  }
  if (!"sample_retained" %in% colnames(thresholds)) thresholds$sample_retained <- TRUE
  safe_write_tsv(thresholds, file.path(output_dir, "sample_specific_qc_thresholds.tsv"))
  basic_state <- list(
    object = object,
    qc_before = qc_before,
    thresholds = thresholds,
    failure_summary = failure_summary
  )
  save_checkpoint_rds(basic_state, basic_checkpoint, step_signature, compress = FALSE)
} else {
  object <- basic_state$object
  qc_before <- basic_state$qc_before
  thresholds <- basic_state$thresholds
  failure_summary <- basic_state$failure_summary
}

doublet_checkpoint <- file.path(checkpoint_dir, "02_after_doublet_filter.rds")
doublet_state <- load_checkpoint_rds(doublet_checkpoint, step_signature)
if (is.null(doublet_state)) {
  doublet_summary <- data.frame()
  object$scDblFinder.class <- "singlet"
  object$scDblFinder.score <- NA_real_

  if (isTRUE(cfg$qc$run_scdblfinder)) {
    require_packages(
      c("SingleCellExperiment", "scDblFinder", "SummarizedExperiment"),
      "sample-level doublet detection"
    )

    all_classes <- setNames(rep(NA_character_, ncol(object)), colnames(object))
    all_scores <- setNames(rep(NA_real_, ncol(object)), colnames(object))
    sample_ids_qc <- sort(unique(as.character(object$sample_id)))

    for (i in seq_along(sample_ids_qc)) {
      sample_id <- sample_ids_qc[[i]]
      cells <- colnames(object)[as.character(object$sample_id) == sample_id]
      log_message(
        "scDblFinder: ", sample_id, " with ", length(cells), " cells (",
        i, "/", length(sample_ids_qc), ")"
      )

      sample_counts <- get_assay_counts(object[, cells, drop = FALSE])
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = sample_counts)
      )
      set.seed(cfg$seed + i)
      sce <- scDblFinder::scDblFinder(sce, verbose = FALSE)
      sce_md <- as.data.frame(SummarizedExperiment::colData(sce))
      all_classes[cells] <- as.character(sce_md$scDblFinder.class)
      all_scores[cells] <- as.numeric(sce_md$scDblFinder.score)

      doublet_summary <- rbind(
        doublet_summary,
        data.frame(
          sample_id = sample_id,
          cells_tested = length(cells),
          singlets = sum(sce_md$scDblFinder.class == "singlet"),
          doublets = sum(sce_md$scDblFinder.class == "doublet"),
          doublet_fraction = mean(sce_md$scDblFinder.class == "doublet"),
          stringsAsFactors = FALSE
        )
      )
      rm(sample_counts, sce, sce_md)
      gc(FALSE)
    }

    object$scDblFinder.class <- all_classes[colnames(object)]
    object$scDblFinder.score <- all_scores[colnames(object)]
    keep_singlets <- object$scDblFinder.class == cfg$qc$keep_doublet_class
    keep_singlets[is.na(keep_singlets)] <- FALSE
    object <- subset(object, cells = colnames(object)[keep_singlets])
  }
  doublet_state <- list(object = object, doublet_summary = doublet_summary)
  save_checkpoint_rds(doublet_state, doublet_checkpoint, step_signature, compress = FALSE)
} else {
  object <- doublet_state$object
  doublet_summary <- doublet_state$doublet_summary
}

safe_write_tsv(doublet_summary, file.path(output_dir, "scDblFinder_by_sample.tsv"))

qc_after <- object@meta.data
qc_after$cell_id <- rownames(qc_after)

summarize_qc <- function(md, stage) {
  result <- stats::aggregate(
    cbind(nFeature_RNA, nCount_RNA, percent.mt) ~ sample_id + subject_id + cohort + condition,
    data = md,
    FUN = function(x) stats::median(x, na.rm = TRUE)
  )
  counts <- as.data.frame(table(md$sample_id), stringsAsFactors = FALSE)
  colnames(counts) <- c("sample_id", "cells")
  result <- merge(result, counts, by = "sample_id", all.x = TRUE)
  result$stage <- stage
  result
}

qc_summary <- rbind(
  summarize_qc(qc_before, "before_qc"),
  summarize_qc(qc_after, "after_qc_doublet")
)
qc_summary <- merge(qc_summary, failure_summary, by = "sample_id", all.x = TRUE)
safe_write_tsv(qc_summary, file.path(output_dir, "qc_summary_by_sample.tsv"))

plot_qc_violin <- function(metadata, title, output_path) {
  plot_data <- data.table::melt(
    data.table::as.data.table(metadata),
    id.vars = "cohort_condition",
    measure.vars = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    variable.name = "metric",
    value.name = "value"
  )
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = cohort_condition, y = value, fill = cohort_condition)
  ) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL)
  ggplot2::ggsave(output_path, plot = plot, width = 14, height = 5, units = "in")
}

plot_qc_violin(
  qc_before,
  "GSE135779 QC before filtering",
  file.path(output_dir, "qc_violin_before.pdf")
)
plot_qc_violin(
  qc_after,
  "GSE135779 QC after filtering and doublet removal",
  file.path(output_dir, "qc_violin_after.pdf")
)

embedding_checkpoint <- file.path(checkpoint_dir, "03_normalized_pca_harmony.rds")
embedding_state <- load_checkpoint_rds(embedding_checkpoint, step_signature)
if (is.null(embedding_state)) {
  log_message("Normalizing data and selecting variable features.")
  object <- Seurat::NormalizeData(
    object,
    normalization.method = "LogNormalize",
    scale.factor = cfg$dimension_reduction$normalization_scale_factor,
    verbose = FALSE
  )
  object <- Seurat::FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = cfg$dimension_reduction$n_variable_features,
    verbose = FALSE
  )
  object <- Seurat::ScaleData(
    object,
    features = Seurat::VariableFeatures(object),
    verbose = FALSE
  )
  object <- Seurat::RunPCA(
    object,
    features = Seurat::VariableFeatures(object),
    npcs = cfg$dimension_reduction$n_pcs,
    seed.use = cfg$seed,
    verbose = FALSE
  )

  downstream_reduction <- "pca"
  if (isTRUE(cfg$dimension_reduction$run_harmony)) {
    require_packages("harmony", "Harmony embedding")
    harmony_group <- cfg$dimension_reduction$harmony_group
    if (!harmony_group %in% colnames(object@meta.data)) {
      stop("Harmony group is absent from cell metadata: ", harmony_group)
    }
    harmony_values <- as.character(object@meta.data[[harmony_group]])
    valid_harmony <- !is.na(harmony_values) & nzchar(harmony_values)
    if (!all(valid_harmony)) stop("Harmony group contains missing values: ", harmony_group)
    if (length(unique(harmony_values)) >= 2L) {
      log_message("Running Harmony on technical field ", harmony_group)
      object <- harmony::RunHarmony(
        object = object,
        group.by.vars = harmony_group,
        reduction.use = "pca",
        dims.use = cfg$dimension_reduction$dims_use,
        plot_convergence = FALSE,
        verbose = FALSE
      )
      downstream_reduction <- "harmony"
    }
  }
  embedding_state <- list(object = object, downstream_reduction = downstream_reduction)
  save_checkpoint_rds(
    embedding_state,
    embedding_checkpoint,
    step_signature,
    compress = FALSE
  )
} else {
  object <- embedding_state$object
  downstream_reduction <- embedding_state$downstream_reduction
}

object <- Seurat::FindNeighbors(
  object,
  reduction = downstream_reduction,
  dims = cfg$dimension_reduction$dims_use,
  verbose = FALSE
)

resolution_summary <- list()
for (resolution in cfg$dimension_reduction$resolution_grid) {
  object <- Seurat::FindClusters(
    object,
    resolution = resolution,
    random.seed = cfg$seed,
    verbose = FALSE
  )
  resolution_summary[[as.character(resolution)]] <- data.frame(
    resolution = resolution,
    n_clusters = length(unique(as.character(SeuratObject::Idents(object)))),
    selected = isTRUE(all.equal(resolution, cfg$dimension_reduction$resolution)),
    stringsAsFactors = FALSE
  )
}
resolution_summary <- do.call(rbind, resolution_summary)
safe_write_tsv(
  resolution_summary,
  file.path(output_dir, "cluster_resolution_summary.tsv")
)
object <- Seurat::FindClusters(
  object,
  resolution = cfg$dimension_reduction$resolution,
  random.seed = cfg$seed,
  verbose = FALSE
)

object <- Seurat::RunUMAP(
  object,
  reduction = downstream_reduction,
  dims = cfg$dimension_reduction$dims_use,
  n.neighbors = cfg$dimension_reduction$umap_neighbors,
  min.dist = cfg$dimension_reduction$umap_min_dist,
  seed.use = cfg$seed,
  verbose = FALSE
)

umap_cluster <- Seurat::DimPlot(
  object,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) + ggplot2::ggtitle("GSE135779 clusters")

umap_condition <- Seurat::DimPlot(
  object,
  reduction = "umap",
  group.by = "condition",
  split.by = "cohort",
  raster = TRUE
) + ggplot2::ggtitle("Condition by cohort")

ggplot2::ggsave(
  file.path(output_dir, "umap_clusters.pdf"),
  umap_cluster,
  width = 9,
  height = 7
)
ggplot2::ggsave(
  file.path(output_dir, "umap_condition_cohort.pdf"),
  umap_condition,
  width = 13,
  height = 6
)

cluster_counts <- as.data.frame(table(
  cluster = object$seurat_clusters,
  cohort = object$cohort,
  condition = object$condition
), stringsAsFactors = FALSE)
safe_write_tsv(cluster_counts, file.path(output_dir, "cluster_cell_counts.tsv"))

object@misc$analysis_config <- cfg
object@misc$downstream_reduction <- downstream_reduction
atomic_save_rds(object, output_file, compress = FALSE)
write_session_info(file.path(cfg$paths$logs, "03_qc_cluster_sessionInfo.txt"))
critical_outputs <- c(
  output_file,
  file.path(output_dir, "sample_specific_qc_thresholds.tsv"),
  file.path(output_dir, "scDblFinder_by_sample.tsv"),
  file.path(output_dir, "qc_summary_by_sample.tsv"),
  file.path(output_dir, "cluster_resolution_summary.tsv"),
  file.path(output_dir, "cluster_cell_counts.tsv"),
  file.path(output_dir, "qc_violin_before.pdf"),
  file.path(output_dir, "qc_violin_after.pdf"),
  file.path(output_dir, "umap_clusters.pdf"),
  file.path(output_dir, "umap_condition_cohort.pdf")
)
complete_step_and_cleanup(
  output_dir,
  step_signature,
  critical_outputs,
  cleanup_paths = checkpoint_dir
)

log_message("QC and clustering object saved: ", output_file)
log_message("Retained ", ncol(object), " cells across ", length(unique(object$sample_id)), " samples.")
