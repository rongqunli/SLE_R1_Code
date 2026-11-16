if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
options(stringsAsFactors = FALSE)

shared_root <- normalizePath(
  file.path(PIPELINE_ROOT, "..", "GSE88884_pipeline"),
  winslash = "/",
  mustWork = TRUE
)
.libPaths(c(file.path(shared_root, "R_library"), .libPaths()))

required_packages <- c(
  "data.table", "matrixStats", "AnnotationDbi",
  "hta20transcriptcluster.db", "ggplot2", "lme4"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(AnnotationDbi)
  library(hta20transcriptcluster.db)
  library(ggplot2)
  library(lme4)
})

downloads_dir <- file.path(PIPELINE_ROOT, "data", "raw", "downloads")
processed_dir <- file.path(PIPELINE_ROOT, "data", "processed")
results_dir <- file.path(PIPELINE_ROOT, "results")
audit_dir <- file.path(results_dir, "00_input_audit")
expression_dir <- file.path(results_dir, "01_target_expression")
score_dir <- file.path(results_dir, "02_program_scores")
model_dir <- file.path(results_dir, "03_models")
summary_dir <- file.path(results_dir, "04_summary")
log_dir <- file.path(results_dir, "logs")
invisible(lapply(
  c(processed_dir, audit_dir, expression_dir, score_dir, model_dir, summary_dir, log_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

series_file <- file.path(
  downloads_dir,
  paste0(cfg$accession, "_series_matrix.txt.gz")
)
expression_file <- file.path(downloads_dir, cfg$expression_file)
contract_file <- file.path(PIPELINE_ROOT, "config", "analysis_contract.tsv")
input_files <- c(series_file, expression_file, contract_file)
if (any(!file.exists(input_files))) {
  stop("Missing required inputs: ", paste(input_files[!file.exists(input_files)], collapse = ", "))
}

strip_geo_quotes <- function(x) {
  gsub('^"|"$', "", x)
}

read_geo_header <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  lines <- character()
  repeat {
    chunk <- readLines(con, n = 256L, warn = FALSE)
    if (!length(chunk)) break
    stop_at <- match("!series_matrix_table_begin", chunk, nomatch = 0L)
    if (stop_at > 0L) {
      if (stop_at > 1L) lines <- c(lines, chunk[seq_len(stop_at - 1L)])
      break
    }
    lines <- c(lines, chunk)
  }
  lines
}

geo_values <- function(header_lines, key) {
  prefix <- paste0(key, "\t")
  hit <- header_lines[startsWith(header_lines, prefix)]
  if (!length(hit)) return(character())
  strip_geo_quotes(strsplit(hit[[1L]], "\t", fixed = TRUE)[[1L]][-1L])
}

header_lines <- read_geo_header(series_file)
sample_id <- geo_values(header_lines, "!Sample_geo_accession")
sample_title <- geo_values(header_lines, "!Sample_title")
source_name <- geo_values(header_lines, "!Sample_source_name_ch1")
supplementary <- geo_values(header_lines, "!Sample_supplementary_file")
if (!identical(length(sample_id), cfg$expected_samples) ||
    any(c(length(sample_title), length(source_name), length(supplementary)) != length(sample_id))) {
  stop("Unexpected GEO sample metadata dimensions.")
}

metadata <- data.table(
  sample_id = sample_id,
  sample_title = sample_title,
  source_name = source_name,
  supplementary_file = supplementary
)
metadata[, subject_id := trimws(sub(",.*$", "", sample_title))]
metadata[, condition := fifelse(
  grepl("normal", source_name, ignore.case = TRUE),
  "Healthy",
  fifelse(grepl("SLE", source_name, ignore.case = TRUE), "SLE", "Unknown")
)]
metadata[, time := tolower(fcase(
  grepl("week16", sample_title, ignore.case = TRUE), "week16",
  grepl("week52", sample_title, ignore.case = TRUE), "week52",
  grepl("baseline", sample_title, ignore.case = TRUE), "baseline",
  default = "unknown"
))]
metadata[, treatment_observed := fcase(
  grepl("Q2W", source_name, ignore.case = TRUE), "Tabalumab_Q2W",
  grepl("Q4W", source_name, ignore.case = TRUE), "Tabalumab_Q4W",
  grepl("Placebo", source_name, ignore.case = TRUE), "Placebo",
  default = "Unknown"
)]
metadata[, treatment := {
  known <- treatment_observed[treatment_observed != "Unknown"]
  if (length(known)) known[[1L]] else "Unknown"
}, by = subject_id]
metadata[, expression_column := basename(sub("\\.gz$", "", supplementary_file))]
metadata[, expression_column := sub("^GSM[0-9]+_", "", expression_column)]

con <- gzfile(expression_file, "rt")
expression_header <- strip_geo_quotes(
  strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]][-1L]
)
close(con)
if (length(expression_header) != cfg$expected_samples ||
    !setequal(expression_header, metadata$expression_column)) {
  stop("Expression columns do not map one-to-one to GEO samples.")
}
metadata <- metadata[match(expression_header, expression_column)]
if (anyNA(metadata$sample_id)) stop("Failed to order metadata to the expression matrix.")

primary_path <- file.path(
  PIPELINE_ROOT, "..", "GSE200306_pipeline", "results",
  "07_cross_dataset_bridge", "expanded_program_membership.tsv"
)
secondary_path <- file.path(
  PIPELINE_ROOT, "..", "GSE65391_pipeline", "config", "programs.tsv"
)
if (!file.exists(primary_path) || !file.exists(secondary_path)) {
  stop("Frozen program definitions are missing.")
}
primary <- fread(primary_path, sep = "\t")
primary <- primary[
  program %chin% c(
    "CYAN_NEUTROPHIL_IMMATURE",
    "GREENYELLOW_B_CELL_COMPLEMENT",
    "MAGENTA_IFN_ANTIVIRAL"
  ),
  .(program, gene = toupper(gene), analysis_role = "primary")
]
secondary <- fread(secondary_path, sep = "\t")
secondary <- secondary[
  program %chin% c("CYTOTOXICITY_SHARED", "CD8_GZMK_STATE"),
  .(program, genes, analysis_role = "secondary")
]
secondary <- secondary[, .(
  gene = toupper(trimws(unlist(strsplit(genes, ";", fixed = TRUE))))
), by = .(program, analysis_role)]
programs <- unique(rbindlist(list(primary, secondary), use.names = TRUE))
programs <- programs[nzchar(gene)]

annotation <- suppressMessages(AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = unique(programs$gene),
  keytype = "SYMBOL",
  columns = "PROBEID"
))
setDT(annotation)
setnames(annotation, c("SYMBOL", "PROBEID"), c("gene", "probe_id"))
annotation[, gene := toupper(gene)]
annotation <- annotation[
  !is.na(probe_id) & startsWith(probe_id, "TC"),
  .(gene, probe_id)
]
annotation <- unique(annotation)
target_probes <- unique(annotation$probe_id)
if (!length(target_probes)) stop("No target transcript clusters were annotated.")

message(cfg$accession, ": streaming target rows from the deposited expression matrix.")
con <- gzfile(expression_file, "rt")
header_line <- readLines(con, n = 1L, warn = FALSE)
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
    if (!probe_id %chin% target_probes) next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    values <- suppressWarnings(as.numeric(strip_geo_quotes(fields[-1L])))
    if (length(values) != cfg$expected_samples || any(!is.finite(values))) {
      stop("Invalid target expression row: ", raw_id)
    }
    matched_values[[length(matched_values) + 1L]] <- values
    matched_probe <- c(matched_probe, probe_id)
    matched_raw_id <- c(matched_raw_id, raw_id)
  }
}
close(con)
if (row_count != cfg$expected_rows) {
  stop("Unexpected expression row count: ", row_count)
}
if (!length(matched_values)) stop("No frozen-program expression rows were found.")

probe_expression <- do.call(rbind, matched_values)
rownames(probe_expression) <- matched_raw_id
colnames(probe_expression) <- metadata$sample_id
probe_audit <- data.table(
  row_index = seq_along(matched_probe),
  probe_id = matched_probe,
  raw_row_id = matched_raw_id,
  variance = matrixStats::rowVars(probe_expression)
)
probe_audit <- merge(probe_audit, annotation, by = "probe_id", allow.cartesian = TRUE)
setorder(probe_audit, gene, -variance, probe_id)
probe_audit[, selected := seq_len(.N) == 1L, by = gene]
selected <- probe_audit[selected == TRUE]
gene_expression <- probe_expression[selected$row_index, , drop = FALSE]
rownames(gene_expression) <- selected$gene

reference_index <- if (identical(cfg$role, "longitudinal")) {
  which(metadata$time == "baseline")
} else {
  seq_len(nrow(metadata))
}
gene_center <- rowMeans(gene_expression[, reference_index, drop = FALSE])
gene_scale <- matrixStats::rowSds(gene_expression[, reference_index, drop = FALSE])
gene_scale[!is.finite(gene_scale) | gene_scale == 0] <- NA_real_
gene_z <- sweep(gene_expression, 1L, gene_center, "-")
gene_z <- sweep(gene_z, 1L, gene_scale, "/")

coverage <- programs[, .(
  n_defined = uniqueN(gene),
  n_matched = uniqueN(gene[gene %chin% rownames(gene_z)]),
  matched_genes = paste(sort(unique(gene[gene %chin% rownames(gene_z)])), collapse = ";")
), by = .(program, analysis_role)]
coverage[, coverage_fraction := n_matched / n_defined]
coverage[, score_eligible := n_matched >= 5L]

score_matrix <- sapply(coverage$program, function(program_name) {
  genes <- intersect(
    programs[program == program_name, unique(gene)],
    rownames(gene_z)
  )
  if (length(genes) < 5L) return(rep(NA_real_, ncol(gene_z)))
  colMeans(gene_z[genes, , drop = FALSE], na.rm = TRUE)
})
if (is.null(dim(score_matrix))) {
  score_matrix <- matrix(score_matrix, ncol = 1L)
  colnames(score_matrix) <- coverage$program[[1L]]
}
scores <- cbind(copy(metadata), as.data.table(score_matrix))
score_long <- melt(
  scores,
  id.vars = names(metadata),
  measure.vars = coverage$program,
  variable.name = "program",
  value.name = "score_z"
)
score_long <- merge(
  score_long,
  coverage[, .(program, analysis_role)],
  by = "program",
  sort = FALSE
)

file_info <- file.info(c(series_file, expression_file))
download_manifest <- data.table(
  file = basename(rownames(file_info)),
  bytes = as.numeric(file_info$size),
  modified_time = as.character(file_info$mtime),
  source = c("NCBI GEO series matrix", "NCBI GEO supplementary processed matrix")
)
fwrite(download_manifest, file.path(audit_dir, "download_manifest.tsv"), sep = "\t")
fwrite(metadata, file.path(audit_dir, "sample_audit.tsv"), sep = "\t")
input_summary <- data.table(
  metric = c(
    "samples", "subjects", "expression_rows", "expression_columns",
    "mapped_expression_columns", "target_transcript_rows", "selected_target_genes"
  ),
  value = c(
    nrow(metadata), uniqueN(metadata$subject_id), row_count, length(expression_header),
    sum(expression_header %chin% metadata$expression_column),
    nrow(probe_expression), nrow(gene_expression)
  )
)
fwrite(input_summary, file.path(audit_dir, "input_audit_summary.tsv"), sep = "\t")
fwrite(probe_audit, file.path(expression_dir, "target_probe_mapping.tsv"), sep = "\t")
fwrite(
  data.table(
    metric = c("target_rows", "selected_genes", "samples", "reference_samples"),
    value = c(nrow(probe_expression), nrow(gene_expression), ncol(gene_expression), length(reference_index))
  ),
  file.path(expression_dir, "expression_dimension_audit.tsv"),
  sep = "\t"
)
fwrite(coverage, file.path(score_dir, "program_coverage.tsv"), sep = "\t")
fwrite(score_long, file.path(score_dir, "sample_program_scores.tsv"), sep = "\t")

extract_coefficient <- function(fit, term, feature, family, n_subjects, n_samples) {
  coefs <- coef(summary(fit))
  result <- data.table(
    model_family = family,
    program = feature,
    term = term,
    estimate = NA_real_,
    std_error = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    n_subjects = n_subjects,
    n_samples = n_samples,
    p_method = "two-sided normal approximation"
  )
  if (!term %in% rownames(coefs)) return(result)
  beta_value <- unname(coefs[term, "Estimate"])
  se_value <- unname(coefs[term, "Std. Error"])
  z_value <- beta_value / se_value
  p_value_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
  result[, `:=`(
    estimate = beta_value,
    std_error = se_value,
    conf_low = beta_value - 1.96 * se_value,
    conf_high = beta_value + 1.96 * se_value,
    statistic = z_value,
    p_value = p_value_value
  )]
  result
}

if (identical(cfg$role, "baseline")) {
  model_results <- rbindlist(lapply(unique(score_long$program), function(program_name) {
    dat <- score_long[program == program_name & is.finite(score_z)]
    dat[, condition := relevel(factor(condition), ref = "Healthy")]
    fit <- lm(score_z ~ condition, data = dat)
    extract_coefficient(
      fit,
      "conditionSLE",
      program_name,
      "baseline_case_control",
      uniqueN(dat$subject_id),
      nrow(dat)
    )
  }))
  model_results <- merge(
    model_results,
    coverage[, .(program, analysis_role)],
    by = "program",
    sort = FALSE
  )
  model_results[, fdr := p.adjust(p_value, method = "BH"), by = analysis_role]
  fwrite(model_results, file.path(model_dir, "baseline_model_results.tsv"), sep = "\t")

  p_scores <- ggplot(
    score_long,
    aes(x = condition, y = score_z, fill = condition)
  ) +
    geom_boxplot(outlier.shape = NA, width = 0.65) +
    facet_wrap(~program, scales = "free_y", ncol = 2L) +
    coord_cartesian(ylim = quantile(score_long$score_z, c(0.01, 0.99), na.rm = TRUE)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none") +
    labs(
      title = paste0(cfg$accession, " frozen program scores"),
      x = NULL,
      y = "Program score (gene-wise z mean)"
    )
  ggsave(
    file.path(summary_dir, "program_case_control.png"),
    p_scores,
    width = 9,
    height = 7,
    dpi = 180
  )
  sensitivity_results <- data.table()
} else {
  score_long[, time := relevel(factor(time), ref = "baseline")]
  score_long[, treatment := relevel(factor(treatment), ref = "Placebo")]
  target_terms <- c(
    "timeweek16",
    "timeweek52",
    "timeweek16:treatmentTabalumab_Q2W",
    "timeweek52:treatmentTabalumab_Q2W",
    "timeweek16:treatmentTabalumab_Q4W",
    "timeweek52:treatmentTabalumab_Q4W"
  )
  model_results <- rbindlist(lapply(unique(score_long$program), function(program_name) {
    dat <- score_long[program == program_name & is.finite(score_z)]
    fit <- lmer(
      score_z ~ time * treatment + (1 | subject_id),
      data = dat,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 200000),
        check.conv.singular = "ignore"
      )
    )
    rbindlist(lapply(target_terms, function(term) {
      extract_coefficient(
        fit,
        term,
        program_name,
        "longitudinal_mixed",
        uniqueN(dat$subject_id),
        nrow(dat)
      )
    }))
  }))
  model_results <- merge(
    model_results,
    coverage[, .(program, analysis_role)],
    by = "program",
    sort = FALSE
  )
  model_results[, fdr := p.adjust(p_value, method = "BH"), by = .(term, analysis_role)]
  fwrite(model_results, file.path(model_dir, "mixed_model_results.tsv"), sep = "\t")

  sensitivity_results <- rbindlist(lapply(unique(score_long$program), function(program_name) {
    dat <- score_long[program == program_name & is.finite(score_z)]
    baseline <- dat[time == "baseline", .(subject_id, baseline_score = score_z)]
    rbindlist(lapply(c("week16", "week52"), function(target_time) {
      followup <- dat[time == target_time, .(subject_id, treatment, followup_score = score_z)]
      paired <- merge(followup, baseline, by = "subject_id")
      paired[, delta := followup_score - baseline_score]
      paired[, treatment := relevel(factor(treatment), ref = "Placebo")]
      fit <- lm(delta ~ treatment, data = paired)
      rbindlist(lapply(c("treatmentTabalumab_Q2W", "treatmentTabalumab_Q4W"), function(term) {
        extract_coefficient(
          fit,
          term,
          program_name,
          paste0("paired_change_", target_time),
          uniqueN(paired$subject_id),
          nrow(paired)
        )
      }))
    }))
  }))
  sensitivity_results <- merge(
    sensitivity_results,
    coverage[, .(program, analysis_role)],
    by = "program",
    sort = FALSE
  )
  sensitivity_results[, fdr := p.adjust(p_value, method = "BH"),
                      by = .(model_family, analysis_role)]
  fwrite(
    sensitivity_results,
    file.path(model_dir, "paired_change_sensitivity.tsv"),
    sep = "\t"
  )

  trajectory <- score_long[, .(
    mean = mean(score_z, na.rm = TRUE),
    se = sd(score_z, na.rm = TRUE) / sqrt(.N),
    n = .N
  ), by = .(program, treatment, time)]
  trajectory[, time := factor(time, levels = c("baseline", "week16", "week52"))]
  p_trajectory <- ggplot(
    trajectory,
    aes(x = time, y = mean, color = treatment, group = treatment)
  ) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0.08) +
    facet_wrap(~program, scales = "free_y", ncol = 2L) +
    theme_bw(base_size = 9) +
    labs(
      title = paste0(cfg$accession, " longitudinal pharmacodynamic scores"),
      x = NULL,
      y = "Mean program score",
      color = "Treatment"
    )
  ggsave(
    file.path(summary_dir, "program_trajectories.png"),
    p_trajectory,
    width = 10,
    height = 7,
    dpi = 180
  )
}

saveRDS(
  list(
    metadata = metadata,
    gene_expression = gene_expression,
    program_scores = score_long,
    coverage = coverage,
    models = model_results,
    sensitivity = sensitivity_results,
    annotation_package = "hta20transcriptcluster.db 8.8.0"
  ),
  file.path(processed_dir, "analysis_results.rds")
)

if (identical(cfg$role, "baseline")) {
  key_table <- model_results[order(fdr), .(
    program, estimate, conf_low, conf_high, p_value, fdr
  )]
  key_lines <- paste0(
    "- ", key_table$program,
    ": SLE−Healthy = ", sprintf("%.3f", key_table$estimate),
    " (95% CI ", sprintf("%.3f", key_table$conf_low), " to ",
    sprintf("%.3f", key_table$conf_high), "), FDR = ",
    format(key_table$fdr, digits = 3)
  )
  interpretation <- c(
    "本数据集只支持基线病例—对照关联分析。",
    "公开数据没有个体级 flare 时间，因此不复现 flare 预测模型。"
  )
} else {
  key_table <- model_results[
    grepl(":", term) & program == "GREENYELLOW_B_CELL_COMPLEMENT"
  ][order(term)]
  key_lines <- paste0(
    "- ", key_table$term,
    ": differential change = ", sprintf("%.3f", key_table$estimate),
    ", FDR = ", format(key_table$fdr, digits = 3)
  )
  interpretation <- c(
    "本数据集检验的是相对安慰剂的纵向药效学改变，不是临床疗效预测。",
    "公开数据没有 responder、SLEDAI 或 flare 终点，因此预测模型按契约停止。"
  )
}
summary_lines <- c(
  paste0("# ", cfg$accession, " 首轮分析结果"),
  "",
  paste0("- 样本数：", nrow(metadata)),
  paste0("- 受试者数：", uniqueN(metadata$subject_id)),
  paste0("- 表达矩阵行数：", row_count),
  paste0("- 冻结程序覆盖的唯一基因数：", nrow(gene_expression)),
  "",
  "## 主要结果",
  "",
  key_lines,
  "",
  "## 解释边界",
  "",
  paste0("- ", interpretation),
  "",
  "所有程序在查看结果前已由既有管线冻结；本管线没有在当前数据集中重新选择基因。"
)
writeLines(
  summary_lines,
  file.path(PIPELINE_ROOT, paste0(cfg$accession, "_result_summary.md")),
  useBytes = TRUE
)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
saveRDS(
  list(
    accession = cfg$accession,
    completed = Sys.time(),
    samples = nrow(metadata),
    subjects = uniqueN(metadata$subject_id)
  ),
  file.path(summary_dir, "_SUCCESS.rds")
)
message(cfg$accession, " analysis complete.")
