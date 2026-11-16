if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
.libPaths(unique(c(
  "F:/program/Hei/SLE/.reviewer_Rlib",
  "F:/program/Hei/SLE_v1/pipeline/3.2/data/annotation_library",
  .libPaths()
)))

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(AnnotationDbi)
  library(SeuratObject)
  library(Matrix)
  library(quantiseqr)
})

project_root <- "F:/program/Hei/SLE_v1"
source_root <- "F:/program/Hei/SLE"
output_root <- file.path(project_root, "results", "reviewer_response", "07_bcell_abundance")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "checkpoints"), recursive = TRUE, showWarnings = FALSE)

set.seed(20260811L)

write_tsv <- function(x, name) {
  fwrite(as.data.table(x), file.path(output_root, name), sep = "\t", na = "NA")
}

strip_geo_quotes <- function(x) {
  x <- sub('^"', "", x)
  sub('"$', "", x)
}

safe_z <- function(x, reference) {
  center <- mean(reference, na.rm = TRUE)
  spread <- sd(reference, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(NA_real_, length(x)))
  (x - center) / spread
}

extract_lm_term <- function(fit, term, label) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) stop("Missing model term: ", term)
  estimate <- unname(sm[term, "Estimate"])
  standard_error <- unname(sm[term, "Std. Error"])
  degrees <- df.residual(fit)
  critical <- qt(0.975, degrees)
  data.table(
    model = label,
    estimate = estimate,
    standard_error = standard_error,
    conf_low = estimate - critical * standard_error,
    conf_high = estimate + critical * standard_error,
    statistic = estimate / standard_error,
    p_value = unname(sm[term, "Pr(>|t|)"])
  )
}

program_membership <- fread(
  file.path(project_root, "pipeline", "3.2", "data", "reference", "expanded_program_membership.tsv"),
  sep = "\t"
)
program_membership[, gene := toupper(gene)]
primary_program_genes <- unique(program_membership$gene)
b_program_genes <- unique(program_membership[
  program == "GREENYELLOW_B_CELL_COMPLEMENT", gene
])







marker_checkpoint <- file.path(output_root, "checkpoints", "healthy_sc_b_specificity.rds")
if (file.exists(marker_checkpoint)) {
  marker_info <- readRDS(marker_checkpoint)
} else {
  seurat_path <- file.path(
    source_root, "GSE135779_pipeline", "data", "processed", "03_program_scored_seurat.rds"
  )
  seurat_object <- readRDS(seurat_path)
  metadata <- seurat_object[[]]
  cell_type <- as.character(metadata$cell_type_analysis)
  broad_type <- rep(NA_character_, nrow(metadata))
  broad_type[cell_type %in% c("B_Naive", "B", "Plasmablast")] <- "B"
  broad_type[cell_type %in% c(
    "CD4_T", "CD4_T_Naive", "CD8_T_Cytotoxic", "CD8_T_GZMK", "CD8_T_Naive"
  )] <- "T"
  broad_type[cell_type == "NK"] <- "NK"
  broad_type[cell_type %in% c("CD14_Monocyte", "CD16_Monocyte", "Monocyte")] <- "Monocyte"
  broad_type[cell_type == "pDC"] <- "pDC"
  broad_type[cell_type == "Platelet"] <- "Platelet"
  broad_type[cell_type == "Erythroid"] <- "Erythroid"

  eligible <- which(metadata$condition == "HD" & !is.na(broad_type))
  balanced_cells <- unlist(lapply(split(eligible, broad_type[eligible]), function(index) {
    sample(index, min(5000L, length(index)), replace = FALSE)
  }), use.names = FALSE)
  balanced_type <- broad_type[balanced_cells]
  normalized <- LayerData(seurat_object[["RNA"]], layer = "data")[, balanced_cells, drop = FALSE]
  type_levels <- sort(unique(balanced_type))
  mean_expression <- sapply(type_levels, function(type) {
    Matrix::rowMeans(normalized[, balanced_type == type, drop = FALSE])
  })
  detection_fraction <- sapply(type_levels, function(type) {
    Matrix::rowMeans(normalized[, balanced_type == type, drop = FALSE] > 0)
  })
  b_column <- match("B", colnames(mean_expression))
  marker_info <- data.table(
    gene = toupper(rownames(mean_expression)),
    mean_B = mean_expression[, b_column],
    pct_B = detection_fraction[, b_column],
    max_mean_nonB = apply(mean_expression[, -b_column, drop = FALSE], 1, max),
    max_pct_nonB = apply(detection_fraction[, -b_column, drop = FALSE], 1, max)
  )
  marker_info[, `:=`(
    b_specificity = mean_B - max_mean_nonB,
    b_detection_ratio = pct_B / pmax(max_pct_nonB, 0.01)
  )]
  saveRDS(marker_info, marker_checkpoint, compress = TRUE)
  rm(seurat_object, normalized, mean_expression, detection_fraction)
  invisible(gc())
}

excluded_state_genes <- c(
  "JUN", "FOS", "FOSB", "JUNB", "JUND", "CD69", "CXCR4", "TSC22D3",
  "TNFRSF13B", "TNFRSF13C", "TNFSF13", "TNFSF13B", "NFKB1", "NFKB2", "REL", "RELA"
)
marker_info[, excluded_frozen_program := gene %in% primary_program_genes]
marker_info[, excluded_family :=
  startsWith(gene, "IGH") |
  startsWith(gene, "IGK") |
  startsWith(gene, "IGL") |
  startsWith(gene, "HLA-") |
  startsWith(gene, "RPL") |
  startsWith(gene, "RPS") |
  startsWith(gene, "MT-") |
  gene %in% c("MZB1", "JCHAIN") |
  gene %in% excluded_state_genes
]
marker_candidates <- marker_info[
  !excluded_frozen_program & !excluded_family &
  pct_B >= 0.15 & b_specificity > 0.15 & b_detection_ratio >= 1.5
]
setorder(marker_candidates, -b_specificity, -b_detection_ratio, gene)
strict_marker_30 <- head(marker_candidates$gene, 30L)
strict_marker_15 <- head(marker_candidates$gene, 15L)
marker_info[, marker_rank := match(gene, strict_marker_30)]
marker_info[, selected_marker_30 := gene %in% strict_marker_30]
marker_info[, selected_marker_15 := gene %in% strict_marker_15]
write_tsv(marker_info, "healthy_sc_b_specificity.tsv")
write_tsv(marker_info[selected_marker_30 == TRUE][order(marker_rank)], "orthogonal_b_marker_panel.tsv")







if (!requireNamespace("hta20transcriptcluster.db", quietly = TRUE)) {
  stop("hta20transcriptcluster.db is unavailable")
}
annotation_db <- hta20transcriptcluster.db::hta20transcriptcluster.db
all_probe_ids <- AnnotationDbi::keys(annotation_db, keytype = "PROBEID")
all_probe_ids <- all_probe_ids[startsWith(all_probe_ids, "TC")]
probe_to_gene <- suppressMessages(AnnotationDbi::select(
  annotation_db,
  keys = all_probe_ids,
  keytype = "PROBEID",
  columns = "SYMBOL"
))
probe_to_gene <- as.data.table(probe_to_gene)
probe_to_gene[, SYMBOL := toupper(trimws(SYMBOL))]
probe_to_gene <- probe_to_gene[!is.na(SYMBOL) & nzchar(SYMBOL)]
probe_to_gene[, n_symbols := uniqueN(SYMBOL), by = PROBEID]
probe_to_gene <- unique(probe_to_gene[n_symbols == 1L, .(probe_id = PROBEID, gene = SYMBOL)])
target_transcript_clusters <- unique(probe_to_gene$probe_id)

trial_specs <- list(
  list(
    trial = "ILLUMINATE-1",
    raw_path = file.path(
      source_root, "GSE88885_pipeline", "data", "raw", "downloads",
      "GSE88885_ILLUMINATE1_pd_preprocessed.txt.gz"
    ),
    score_path = file.path(
      project_root, "pipeline", "3.4", "data", "processed",
      "ILLUMINATE1_sample_program_scores.tsv"
    ),
    expected_samples = 908L,
    expected_rows = 92869L
  ),
  list(
    trial = "ILLUMINATE-2",
    raw_path = file.path(
      source_root, "GSE88886_pipeline", "data", "raw", "downloads",
      "GSE88886_ILLUMINATE2_pd_preprocessed.txt.gz"
    ),
    score_path = file.path(
      project_root, "pipeline", "3.4", "data", "processed",
      "ILLUMINATE2_sample_program_scores.tsv"
    ),
    expected_samples = 418L,
    expected_rows = 92869L
  )
)

read_score_metadata <- function(path) {
  scores <- fread(path, sep = "\t")
  sample_metadata <- unique(scores[
    analysis_role == "primary",
    .(sample_id, subject_id, time, treatment, expression_column)
  ])
  b_scores <- scores[
    analysis_role == "primary" & program == "GREENYELLOW_B_CELL_COMPLEMENT",
    .(sample_id, official_b_score = score_z)
  ]
  list(metadata = sample_metadata, b_scores = b_scores)
}

stream_full_gene_expression <- function(spec, sample_metadata) {
  message(spec$trial, ": streaming annotated transcript clusters from the deposited matrix")
  con <- gzfile(spec$raw_path, "rt")
  header_line <- readLines(con, n = 1L, warn = FALSE)
  header_fields <- strip_geo_quotes(strsplit(header_line, "\t", fixed = TRUE)[[1L]])
  expression_columns <- header_fields[-1L]
  if (length(expression_columns) != spec$expected_samples) {
    stop(spec$trial, ": unexpected sample-column count")
  }
  metadata_index <- match(expression_columns, sample_metadata$expression_column)
  if (anyNA(metadata_index)) {
    stop(spec$trial, ": expression columns were not mapped to sample metadata")
  }
  ordered_metadata <- sample_metadata[metadata_index]

  matched_values <- list()
  matched_probe <- character()
  matched_raw_id <- character()
  row_count <- 0L
  repeat {
    chunk <- readLines(con, n = 64L, warn = FALSE)
    if (!length(chunk)) break
    for (line in chunk) {
      row_count <- row_count + 1L
      first_tab <- regexpr("\t", line, fixed = TRUE)[[1L]]
      if (first_tab < 1L) next
      raw_id <- strip_geo_quotes(substr(line, 1L, first_tab - 1L))
      probe_id <- sub(":[0-9]+$", "", raw_id)
      if (!probe_id %chin% target_transcript_clusters) next
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      values <- suppressWarnings(as.numeric(strip_geo_quotes(fields[-1L])))
      if (length(values) != spec$expected_samples || any(!is.finite(values))) {
        stop(spec$trial, ": invalid expression row: ", raw_id)
      }
      matched_values[[length(matched_values) + 1L]] <- values
      matched_probe <- c(matched_probe, probe_id)
      matched_raw_id <- c(matched_raw_id, raw_id)
    }
  }
  close(con)
  if (row_count != spec$expected_rows) {
    stop(spec$trial, ": unexpected expression-row count: ", row_count)
  }

  probe_expression <- do.call(rbind, matched_values)
  rownames(probe_expression) <- matched_raw_id
  colnames(probe_expression) <- ordered_metadata$sample_id
  probe_audit <- data.table(
    row_index = seq_along(matched_probe),
    probe_id = matched_probe,
    raw_row_id = matched_raw_id,
    variance = matrixStats::rowVars(probe_expression)
  )
  probe_audit <- merge(probe_audit, probe_to_gene, by = "probe_id", allow.cartesian = TRUE)
  setorder(probe_audit, gene, -variance, probe_id)
  probe_audit[, selected := seq_len(.N) == 1L, by = gene]
  selected <- probe_audit[selected == TRUE]
  gene_expression <- probe_expression[selected$row_index, , drop = FALSE]
  rownames(gene_expression) <- selected$gene
  if (anyDuplicated(rownames(gene_expression))) stop("Duplicated gene symbols after selection")
  rm(probe_expression, matched_values)
  invisible(gc())
  list(expression = gene_expression, metadata = ordered_metadata, probe_audit = selected)
}

score_gene_set <- function(expression, metadata, genes) {
  genes <- intersect(genes, rownames(expression))
  if (length(genes) < 5L) stop("Fewer than five genes were available for a score")
  baseline <- metadata$time == "baseline"
  centers <- rowMeans(expression[genes, baseline, drop = FALSE])
  scales <- matrixStats::rowSds(expression[genes, baseline, drop = FALSE])
  valid <- is.finite(scales) & scales > 0
  genes <- genes[valid]
  z <- sweep(expression[genes, , drop = FALSE], 1L, centers[valid], "-")
  z <- sweep(z, 1L, scales[valid], "/")
  list(score = colMeans(z, na.rm = TRUE), genes = genes, z = z)
}

sample_tables <- list()
program_gene_effect_tables <- list()
trial_audits <- list()

for (spec in trial_specs) {
  score_input <- read_score_metadata(spec$score_path)
  full <- stream_full_gene_expression(spec, score_input$metadata)

  message(spec$trial, ": running quanTIseq in microarray mode")
  quantiseq_result <- quantiseqr::run_quantiseq(
    full$expression,
    signature_matrix = "TIL10",
    is_arraydata = TRUE,
    is_tumordata = FALSE,
    scale_mRNA = TRUE,
    method = "lsei"
  )
  quantiseq_result <- as.data.table(quantiseq_result)
  b_fraction_column <- grep("^B", names(quantiseq_result), value = TRUE)[1L]
  if (is.na(b_fraction_column)) stop("quanTIseq B-cell column was not found")
  setnames(quantiseq_result, c("Sample", b_fraction_column), c("sample_id", "quantiseq_b_fraction"))
  quantiseq_result <- quantiseq_result[, .(sample_id, quantiseq_b_fraction)]

  proxy30 <- score_gene_set(full$expression, full$metadata, strict_marker_30)
  proxy15 <- score_gene_set(full$expression, full$metadata, strict_marker_15)
  b_program <- score_gene_set(full$expression, full$metadata, b_program_genes)

  sample_table <- copy(full$metadata)
  sample_table[, `:=`(
    marker_proxy_30 = proxy30$score,
    marker_proxy_15 = proxy15$score,
    recomputed_b_program = b_program$score,
    trial = spec$trial
  )]
  sample_table <- merge(sample_table, score_input$b_scores, by = "sample_id", all.x = TRUE, sort = FALSE)
  sample_table <- merge(sample_table, quantiseq_result, by = "sample_id", all.x = TRUE, sort = FALSE)
  baseline_index <- sample_table$time == "baseline"
  sample_table[, quantiseq_b_fraction_z := safe_z(
    quantiseq_b_fraction, quantiseq_b_fraction[baseline_index]
  )]
  sample_table[, marker_proxy_30_z := safe_z(marker_proxy_30, marker_proxy_30[baseline_index])]
  sample_table[, marker_proxy_15_z := safe_z(marker_proxy_15, marker_proxy_15[baseline_index])]
  sample_tables[[spec$trial]] <- sample_table

  sc_program <- marker_info[match(b_program$genes, gene)]
  sc_program[, b_identity_high :=
    is.finite(b_specificity) & b_specificity > 0.10 &
    is.finite(b_detection_ratio) & b_detection_ratio > 1.20
  ]
  high_genes <- sc_program[b_identity_high == TRUE, gene]
  low_genes <- sc_program[b_identity_high == FALSE, gene]
  high_score <- if (length(high_genes) >= 5L) score_gene_set(full$expression, full$metadata, high_genes) else NULL
  low_score <- if (length(low_genes) >= 5L) score_gene_set(full$expression, full$metadata, low_genes) else NULL
  if (!is.null(high_score)) sample_table[, b_identity_high_score := high_score$score]
  if (!is.null(low_score)) sample_table[, b_identity_low_score := low_score$score]
  sample_tables[[spec$trial]] <- sample_table

  validation_cor <- cor(
    sample_table$official_b_score,
    sample_table$recomputed_b_program,
    use = "complete.obs",
    method = "pearson"
  )
  trial_audits[[spec$trial]] <- data.table(
    trial = spec$trial,
    full_genes = nrow(full$expression),
    quantiseq_samples = sum(is.finite(sample_table$quantiseq_b_fraction)),
    marker30_matched = length(proxy30$genes),
    marker15_matched = length(proxy15$genes),
    b_program_matched = length(b_program$genes),
    official_recomputed_b_score_pearson = validation_cor,
    b_identity_high_genes = length(high_genes),
    b_identity_low_genes = length(low_genes)
  )

  
  b_gene_z <- b_program$z
  b_gene_long <- list()
  for (visit in c("week16", "week52")) {
    baseline_cols <- which(full$metadata$time == "baseline")
    followup_cols <- which(full$metadata$time == visit)
    baseline_meta <- full$metadata[baseline_cols]
    followup_meta <- full$metadata[followup_cols]
    baseline_lookup <- match(followup_meta$subject_id, baseline_meta$subject_id)
    complete <- which(!is.na(baseline_lookup))
    delta <- b_gene_z[, followup_cols[complete], drop = FALSE] -
      b_gene_z[, baseline_cols[baseline_lookup[complete]], drop = FALSE]
    delta_meta <- followup_meta[complete]
    for (dose in c("Tabalumab_Q2W", "Tabalumab_Q4W")) {
      drug <- delta_meta$treatment == dose
      placebo <- delta_meta$treatment == "Placebo"
      effect <- rowMeans(delta[, drug, drop = FALSE]) - rowMeans(delta[, placebo, drop = FALSE])
      row <- data.table(
        trial = spec$trial,
        time = visit,
        dose = dose,
        gene = rownames(delta),
        drug_minus_placebo_delta_z = effect
      )
      row <- merge(row, marker_info[, .(gene, b_specificity, b_detection_ratio)], by = "gene", all.x = TRUE)
      b_gene_long[[length(b_gene_long) + 1L]] <- row
    }
  }
  program_gene_effect_tables[[spec$trial]] <- rbindlist(b_gene_long)

  rm(full, quantiseq_result, b_gene_z)
  invisible(gc())
}

sample_level <- rbindlist(sample_tables, fill = TRUE)
write_tsv(sample_level, "sample_level_b_abundance_and_program_scores.tsv")
write_tsv(rbindlist(trial_audits), "analysis_audit.tsv")







correlation_rows <- list()
for (trial_name in unique(sample_level$trial)) {
  for (subset_name in c("baseline", "all_visits")) {
    dat <- if (subset_name == "baseline") {
      sample_level[trial == trial_name & time == "baseline"]
    } else {
      sample_level[trial == trial_name]
    }
    for (proxy in c("quantiseq_b_fraction_z", "marker_proxy_30_z", "marker_proxy_15_z")) {
      for (method in c("spearman", "pearson")) {
        complete <- is.finite(dat$official_b_score) & is.finite(dat[[proxy]])
        test <- cor.test(dat$official_b_score[complete], dat[[proxy]][complete], method = method)
        correlation_rows[[length(correlation_rows) + 1L]] <- data.table(
          trial = trial_name,
          subset = subset_name,
          proxy = proxy,
          method = method,
          n = sum(complete),
          correlation = unname(test$estimate),
          p_value = test$p.value
        )
      }
    }
  }
}
correlation_table <- rbindlist(correlation_rows)
correlation_table[, fdr := p.adjust(p_value, method = "BH")]
write_tsv(correlation_table, "b_program_abundance_correlations.tsv")

make_paired <- function(data, trial_name, visit) {
  baseline <- data[trial == trial_name & time == "baseline"]
  followup <- data[trial == trial_name & time == visit]
  keep_columns <- c(
    "subject_id", "treatment", "official_b_score", "quantiseq_b_fraction_z",
    "marker_proxy_30_z", "marker_proxy_15_z", "b_identity_high_score", "b_identity_low_score"
  )
  keep_columns <- keep_columns[keep_columns %in% names(data)]
  baseline <- baseline[, ..keep_columns]
  followup <- followup[, ..keep_columns]
  setnames(baseline, setdiff(keep_columns, c("subject_id", "treatment")), paste0("baseline_", setdiff(keep_columns, c("subject_id", "treatment"))))
  setnames(followup, setdiff(keep_columns, c("subject_id", "treatment")), paste0("followup_", setdiff(keep_columns, c("subject_id", "treatment"))))
  joined <- merge(baseline, followup, by = c("subject_id", "treatment"), all = FALSE)
  for (column in setdiff(keep_columns, c("subject_id", "treatment"))) {
    joined[, (paste0("delta_", column)) := get(paste0("followup_", column)) - get(paste0("baseline_", column))]
  }
  joined[, `:=`(trial = trial_name, time = visit)]
  joined
}

paired_data <- rbindlist(lapply(unique(sample_level$trial), function(trial_name) {
  rbindlist(lapply(c("week16", "week52"), function(visit) {
    make_paired(sample_level, trial_name, visit)
  }), fill = TRUE)
}), fill = TRUE)
write_tsv(paired_data, "paired_subject_changes.tsv")

bootstrap_coefficient <- function(data, formula, term = "treated", n_boot = 1000L, seed = 1L) {
  set.seed(seed)
  drug <- data[treated == 1L]
  placebo <- data[treated == 0L]
  draws <- rep(NA_real_, n_boot)
  for (iteration in seq_len(n_boot)) {
    boot <- rbind(
      drug[sample(seq_len(nrow(drug)), replace = TRUE)],
      placebo[sample(seq_len(nrow(placebo)), replace = TRUE)]
    )
    fit <- try(lm(formula, data = boot), silent = TRUE)
    if (!inherits(fit, "try-error") && term %in% names(coef(fit))) {
      draws[iteration] <- unname(coef(fit)[term])
    }
  }
  draws <- draws[is.finite(draws)]
  data.table(
    bootstrap_n = length(draws),
    bootstrap_conf_low = quantile(draws, 0.025),
    bootstrap_conf_high = quantile(draws, 0.975),
    bootstrap_p = min(1, 2 * min(mean(draws <= 0), mean(draws >= 0)))
  )
}

model_rows <- list()
seed_counter <- 20260811L
for (trial_name in unique(paired_data$trial)) {
  for (visit in c("week16", "week52")) {
    for (dose in c("Tabalumab_Q2W", "Tabalumab_Q4W")) {
      dat <- copy(paired_data[
        trial == trial_name & time == visit & treatment %in% c("Placebo", dose)
      ])
      dat[, treated := as.integer(treatment == dose)]
      formulas <- list(
        unadjusted_change = delta_official_b_score ~ treated,
        baseline_adjusted_change = delta_official_b_score ~ treated + baseline_official_b_score,
        quantiseq_adjusted_change = delta_official_b_score ~ treated + baseline_official_b_score + baseline_quantiseq_b_fraction_z + delta_quantiseq_b_fraction_z,
        marker30_adjusted_change = delta_official_b_score ~ treated + baseline_official_b_score + baseline_marker_proxy_30_z + delta_marker_proxy_30_z,
        quantiseq_fraction_change = delta_quantiseq_b_fraction_z ~ treated + baseline_quantiseq_b_fraction_z,
        marker30_proxy_change = delta_marker_proxy_30_z ~ treated + baseline_marker_proxy_30_z
      )
      for (model_name in names(formulas)) {
        formula <- formulas[[model_name]]
        needed <- all.vars(formula)
        model_data <- dat[complete.cases(dat[, ..needed])]
        fit <- lm(formula, data = model_data)
        row <- extract_lm_term(fit, "treated", model_name)
        seed_counter <- seed_counter + 1L
        boot <- bootstrap_coefficient(model_data, formula, n_boot = 1000L, seed = seed_counter)
        row <- cbind(row, boot)
        row[, `:=`(
          trial = trial_name,
          time = visit,
          dose = dose,
          n_drug = sum(model_data$treated == 1L),
          n_placebo = sum(model_data$treated == 0L)
        )]
        model_rows[[length(model_rows) + 1L]] <- row
      }
    }
  }
}
model_table <- rbindlist(model_rows, fill = TRUE)
model_table[, fdr := p.adjust(p_value, method = "BH"), by = model]
write_tsv(model_table, "abundance_adjusted_treatment_models.tsv")

effect_wide <- dcast(
  model_table[model %in% c("unadjusted_change", "quantiseq_adjusted_change", "marker30_adjusted_change")],
  trial + time + dose ~ model,
  value.var = "estimate"
)
effect_wide[, quantiseq_attenuation_fraction :=
  1 - quantiseq_adjusted_change / unadjusted_change
]
effect_wide[, marker30_attenuation_fraction :=
  1 - marker30_adjusted_change / unadjusted_change
]
write_tsv(effect_wide, "treatment_effect_attenuation.tsv")


program_gene_effects <- rbindlist(program_gene_effect_tables)
write_tsv(program_gene_effects, "b_program_gene_level_treatment_effects.tsv")
gene_correlation_rows <- program_gene_effects[
  is.finite(b_specificity),
  {
    test <- cor.test(drug_minus_placebo_delta_z, b_specificity, method = "spearman")
    .(
      n_genes = .N,
      rho = unname(test$estimate),
      p_value = test$p.value,
      mean_effect_high_specificity = mean(drug_minus_placebo_delta_z[b_specificity > 0.10]),
      mean_effect_low_specificity = mean(drug_minus_placebo_delta_z[b_specificity <= 0.10])
    )
  },
  by = .(trial, time, dose)
]
gene_correlation_rows[, fdr := p.adjust(p_value, method = "BH")]
write_tsv(gene_correlation_rows, "b_specificity_treatment_effect_correlation.tsv")


baseline_summary <- correlation_table[
  subset == "baseline" & method == "spearman" &
  proxy %in% c("quantiseq_b_fraction_z", "marker_proxy_30_z")
]
primary_models <- model_table[
  model %in% c(
    "unadjusted_change", "quantiseq_adjusted_change", "marker30_adjusted_change",
    "quantiseq_fraction_change"
  )
]
decision_summary <- data.table(
  metric = c(
    "median_baseline_spearman_quantiseq",
    "median_baseline_spearman_marker30",
    "negative_unadjusted_contrasts",
    "negative_quantiseq_adjusted_contrasts",
    "significant_quantiseq_adjusted_contrasts_fdr_0.05",
    "negative_marker30_adjusted_contrasts",
    "significant_marker30_adjusted_contrasts_fdr_0.05",
    "significant_quantiseq_fraction_reductions_fdr_0.05",
    "median_quantiseq_attenuation_fraction",
    "median_marker30_attenuation_fraction"
  ),
  value = c(
    median(baseline_summary[proxy == "quantiseq_b_fraction_z", correlation]),
    median(baseline_summary[proxy == "marker_proxy_30_z", correlation]),
    sum(primary_models[model == "unadjusted_change", estimate < 0]),
    sum(primary_models[model == "quantiseq_adjusted_change", estimate < 0]),
    sum(primary_models[model == "quantiseq_adjusted_change", estimate < 0 & fdr < 0.05]),
    sum(primary_models[model == "marker30_adjusted_change", estimate < 0]),
    sum(primary_models[model == "marker30_adjusted_change", estimate < 0 & fdr < 0.05]),
    sum(primary_models[model == "quantiseq_fraction_change", estimate < 0 & fdr < 0.05]),
    median(effect_wide$quantiseq_attenuation_fraction, na.rm = TRUE),
    median(effect_wide$marker30_attenuation_fraction, na.rm = TRUE)
  )
)
write_tsv(decision_summary, "decision_summary.tsv")

saveRDS(
  list(
    marker_panel_30 = strict_marker_30,
    marker_panel_15 = strict_marker_15,
    correlations = correlation_table,
    treatment_models = model_table,
    attenuation = effect_wide,
    decision_summary = decision_summary
  ),
  file.path(output_root, "reviewer_bcell_abundance_results.rds"),
  compress = TRUE
)

writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
message("Reviewer B-cell abundance analysis completed: ", output_root)
