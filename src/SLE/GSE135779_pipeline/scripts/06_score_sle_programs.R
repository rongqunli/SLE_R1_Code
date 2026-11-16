if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[[1]]), winslash = "/", mustWork = TRUE)), "_common.R"))
analysis_exclusions_file <- file.path(PIPELINE_ROOT, "config", "analysis_exclusions.R")
analysis_exclusion_helpers_file <- file.path(PIPELINE_ROOT, "scripts", "_analysis_exclusion_helpers.R")
source(analysis_exclusions_file, local = FALSE)
source(analysis_exclusion_helpers_file, local = FALSE)

require_packages(
  c("SeuratObject", "Matrix", "data.table", "UCell", "ggplot2"),
  "SLE program scoring"
)

input_file <- file.path(cfg$paths$processed, "02_annotated_seurat.rds")
output_file <- file.path(cfg$paths$processed, "03_program_scored_seurat.rds")
output_dir <- file.path(cfg$paths$results, "06_program_scores")
checkpoint_dir <- file.path(cfg$paths$processed, "checkpoints", "06_program_scores")
score_checkpoint <- file.path(checkpoint_dir, "ucell_score_matrix.rds")
program_file <- file.path(PIPELINE_ROOT, "config", "sle_programs.tsv")
cell_exclusion_file <- file.path(output_dir, "analysis_cell_exclusions.tsv")
cell_exclusion_summary_file <- file.path(output_dir, "analysis_cell_exclusion_summary.tsv")

if (!file.exists(input_file)) {
  stop("Annotated object not found: ", input_file, ". Run step 04_annotate first.")
}

step_signature <- build_step_signature(
  input_files = c(input_file, program_file),
  extra = list(
    step = "06_program_scores",
    programs = cfg$programs,
    pseudobulk = cfg$pseudobulk,
    analysis_exclusions = analysis_exclusions
  ),
  code_files = c(
    script_file(),
    file.path(PIPELINE_ROOT, "scripts", "_common.R"),
    file.path(PIPELINE_ROOT, "config", "config.R"),
    analysis_exclusions_file,
    analysis_exclusion_helpers_file
  )
)
initialize_step(
  step_signature,
  checkpoint_dir,
  owned_paths = c(output_dir, output_file)
)
if (step_is_complete(output_dir, step_signature, c(
  output_file,
  score_checkpoint,
  file.path(output_dir, "program_gene_matching.tsv"),
  file.path(output_dir, "donor_celltype_program_scores.tsv"),
  file.path(output_dir, "donor_level_program_models.tsv"),
  cell_exclusion_file,
  cell_exclusion_summary_file,
  file.path(output_dir, "child_sle_vs_hd_program_score_heatmap.pdf")
))) {
  log_message("Step 06 is already complete and its signature matches. Use --force to rebuild.")
  quit(save = "no", status = 0L)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(input_file)
if (!"cell_type_analysis" %in% colnames(object@meta.data)) {
  stop("cell_type_analysis is missing. Re-run step 04_annotate with the updated annotation code.")
}
analysis_metadata <- object@meta.data
analysis_metadata$cell_id <- rownames(analysis_metadata)
exclusion_result <- mask_excluded_cluster_labels(
  metadata = analysis_metadata,
  excluded_clusters = analysis_exclusions$seurat_clusters,
  label_columns = c("cell_type_analysis"),
  exclusion_reason = analysis_exclusions$reason
)
analysis_metadata <- exclusion_result$metadata
excluded_metadata <- exclusion_result$excluded_metadata

expected_exclusion_counts <- analysis_exclusions$expected_cluster_counts %||% integer(0)
actual_exclusion_counts <- table(as.character(excluded_metadata$seurat_clusters))
if (length(expected_exclusion_counts) > 0L) {
  expected_names <- names(expected_exclusion_counts)
  if (is.null(expected_names) || any(!nzchar(expected_names))) {
    stop("analysis_exclusions$expected_cluster_counts must be a named vector.")
  }
  observed <- as.integer(actual_exclusion_counts[expected_names])
  observed[is.na(observed)] <- 0L
  mismatch <- observed != as.integer(expected_exclusion_counts)
  if (any(mismatch)) {
    stop(
      "Frozen cluster exclusion count mismatch: ",
      paste0(
        expected_names[mismatch], " expected=", as.integer(expected_exclusion_counts[mismatch]),
        " observed=", observed[mismatch], collapse = "; "
      ),
      ". Review the annotated object before running step 06."
    )
  }
}

exclusion_columns <- intersect(
  c(
    "cell_id", "seurat_clusters", "sample_id", "subject_id", "cohort", "condition",
    "cell_type_major", "cell_type_analysis", "analysis_exclusion_reason"
  ),
  colnames(excluded_metadata)
)
safe_write_tsv(excluded_metadata[, exclusion_columns, drop = FALSE], cell_exclusion_file)
exclusion_summary <- data.frame(
  seurat_clusters = names(actual_exclusion_counts),
  n_cells = as.integer(actual_exclusion_counts),
  expected_n_cells = as.integer(expected_exclusion_counts[names(actual_exclusion_counts)]),
  exclusion_reason = analysis_exclusions$reason,
  stringsAsFactors = FALSE
)
safe_write_tsv(exclusion_summary, cell_exclusion_summary_file)
object$analysis_excluded <- exclusion_result$excluded
object$analysis_exclusion_reason <- analysis_metadata$analysis_exclusion_reason
log_message(
  "Excluded ", sum(exclusion_result$excluded), " cells from primary step 06 program models: cluster(s) ",
  paste(analysis_exclusions$seurat_clusters, collapse = ", ")
)
expression_matrix <- get_assay_counts(object)
programs <- read_program_table(program_file)

matching <- do.call(rbind, lapply(names(programs), function(program_name) {
  requested <- unique(programs[[program_name]])
  matched <- intersect(requested, rownames(expression_matrix))
  data.frame(
    program = program_name,
    requested_genes = length(requested),
    matched_genes = length(matched),
    match_fraction = length(matched) / length(requested),
    matched_gene_list = paste(matched, collapse = ";"),
    missing_gene_list = paste(setdiff(requested, matched), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
keep_programs <- matching$program[matching$matched_genes >= cfg$programs$min_matched_genes]
if (length(keep_programs) == 0L) stop("No SLE program has enough matched genes.")
programs <- programs[keep_programs]
programs <- lapply(programs, intersect, y = rownames(expression_matrix))

score_state <- load_checkpoint_rds(score_checkpoint, step_signature)
if (is.null(score_state)) {
  log_message("Calculating rank-based UCell scores for ", length(programs), " programs.")
  score_matrix <- UCell::ScoreSignatures_UCell(expression_matrix, features = programs)
  score_matrix <- as.data.frame(score_matrix, check.names = FALSE)
  if (is.null(rownames(score_matrix)) || !all(colnames(object) %in% rownames(score_matrix))) {
    stop("UCell score rows cannot be matched to Seurat cell names.")
  }
  score_matrix <- score_matrix[colnames(object), , drop = FALSE]
  score_state <- list(score_matrix = score_matrix)
  save_checkpoint_rds(score_state, score_checkpoint, step_signature, compress = TRUE)
} else {
  score_matrix <- score_state$score_matrix
}

original_score_names <- colnames(score_matrix)
program_names_from_ucell <- sub("_UCell$", "", original_score_names)
score_columns <- paste0("program_", make_safe_name(program_names_from_ucell), "_UCell")
if (anyDuplicated(score_columns)) stop("Program names become duplicated after safe-name conversion.")
colnames(score_matrix) <- score_columns
object <- SeuratObject::AddMetaData(object, metadata = score_matrix)
program_column_map <- data.frame(
  program = program_names_from_ucell,
  ucell_output_column = original_score_names,
  seurat_metadata_column = score_columns,
  stringsAsFactors = FALSE
)
matching <- merge(matching, program_column_map, by = "program", all.x = TRUE, sort = FALSE)
safe_write_tsv(matching, file.path(output_dir, "program_gene_matching.tsv"))

score_long <- data.table::melt(
  data.table::as.data.table(data.frame(
    sample_id = analysis_metadata$sample_id,
    subject_id = analysis_metadata$subject_id,
    cohort = analysis_metadata$cohort,
    condition = analysis_metadata$condition,
    age = numeric_or_na(analysis_metadata$age),
    sex = analysis_metadata$sex,
    batch = analysis_metadata$batch,
    SLEDAI = numeric_or_na(analysis_metadata$SLEDAI),
    published_activity_group = analysis_metadata$sle_activity,
    cell_type_analysis = analysis_metadata$cell_type_analysis,
    score_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )),
  id.vars = c(
    "sample_id", "subject_id", "cohort", "condition", "age", "sex", "batch",
    "SLEDAI", "published_activity_group", "cell_type_analysis"
  ),
  measure.vars = score_columns,
  variable.name = "score_column",
  value.name = "score"
)
score_long <- score_long[
  !is.na(cell_type_analysis) & nzchar(cell_type_analysis) &
    !grepl("Unresolved", cell_type_analysis, ignore.case = TRUE)
]
score_long[, program := program_column_map$program[
  match(as.character(score_column), program_column_map$seurat_metadata_column)
]]

donor_scores <- score_long[, .(
  n_cells = .N,
  mean_score = mean(score, na.rm = TRUE),
  median_score = stats::median(score, na.rm = TRUE)
), by = .(
  sample_id, subject_id, cohort, condition, age, sex, batch, SLEDAI,
  published_activity_group, cell_type_analysis, program
)]
donor_scores <- donor_scores[n_cells >= cfg$programs$min_cells_per_donor]
donor_scores[, analysis_activity_group := derive_sledai_group(SLEDAI)]
donor_scores[, activity_binary := data.table::fifelse(
  analysis_activity_group %in% c("Low", "High"), analysis_activity_group, NA_character_
)]
safe_write_tsv(donor_scores, file.path(output_dir, "donor_celltype_program_scores.tsv"))

coerce_covariate <- function(x, name) {
  if (name %in% c("age", "SLEDAI")) return(numeric_or_na(x))
  factor(as.character(x))
}

build_design <- function(
    data,
    effect_column,
    effect_type,
    reference_group = NULL,
    test_group = NULL,
    covariates = character(0)) {
  x <- data
  if (effect_type == "binary") {
    x$.effect <- factor(as.character(x[[effect_column]]), levels = c(reference_group, test_group))
  } else {
    x$.effect <- numeric_or_na(x[[effect_column]])
  }
  x <- x[!is.na(x$.effect), , drop = FALSE]
  if (nrow(x) == 0L) return(NULL)
  used <- character(0)
  final_data <- x
  final_design <- stats::model.matrix(~.effect, data = final_data)
  for (covariate in covariates) {
    if (!covariate %in% colnames(x)) next
    candidate <- x
    candidate_used <- c(used, covariate)
    for (name in candidate_used) candidate[[name]] <- coerce_covariate(candidate[[name]], name)
    complete_columns <- candidate[, c(candidate_used, ".effect"), with = FALSE]
    candidate <- candidate[stats::complete.cases(complete_columns)]
    if (nrow(candidate) == 0L || any(vapply(
      candidate[, candidate_used, with = FALSE],
      function(value) length(unique(value)) < 2L,
      logical(1)
    ))) next
    design <- tryCatch(
      stats::model.matrix(stats::reformulate(c(candidate_used, ".effect")), data = candidate),
      error = function(e) NULL
    )
    if (is.null(design) || qr(design)$rank < ncol(design)) next
    used <- candidate_used
    final_data <- candidate
    final_design <- design
  }
  effect_columns <- grep("^\\.effect", colnames(final_design))
  if (length(effect_columns) != 1L) return(NULL)
  list(
    data = final_data,
    design = final_design,
    coefficient = effect_columns[[1]],
    covariates_used = used,
    formula = paste("~", paste(c(used, effect_column), collapse = " + "))
  )
}

fit_linear <- function(outcome, design, coefficient) {
  inverse <- tryCatch(solve(crossprod(design)), error = function(e) NULL)
  if (is.null(inverse)) return(NULL)
  beta <- as.numeric(inverse %*% crossprod(design, outcome))
  residual <- outcome - as.numeric(design %*% beta)
  df <- nrow(design) - ncol(design)
  if (df <= 0L) return(NULL)
  sigma2 <- sum(residual^2) / df
  se <- sqrt(diag(inverse) * sigma2)
  statistic <- beta[[coefficient]] / se[[coefficient]]
  critical <- stats::qt((1 + cfg$inference$confidence_level) / 2, df = df)
  list(
    estimate = beta[[coefficient]],
    se = se[[coefficient]],
    statistic = statistic,
    df = df,
    p = 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE),
    low = beta[[coefficient]] - critical * se[[coefficient]],
    high = beta[[coefficient]] + critical * se[[coefficient]]
  )
}

run_program_model <- function(
    data,
    comparison,
    effect_column,
    effect_type,
    covariates,
    reference_group = NULL,
    test_group = NULL,
    cohort,
    require_condition = NULL) {
  
  cohort_value <- cohort
  x <- data[as.character(data$cohort) == cohort_value]
  if (!is.null(require_condition)) {
    condition_value <- require_condition
    x <- x[as.character(x$condition) == condition_value]
  }
  results <- x[, {
    design_info <- build_design(
      .SD, effect_column, effect_type, reference_group, test_group, covariates
    )
    if (is.null(design_info)) {
      NULL
    } else {
      if (effect_type == "binary") {
        group_counts <- table(factor(
          design_info$data[[effect_column]], levels = c(reference_group, test_group)
        ))
        enough <- all(group_counts >= cfg$programs$min_samples_per_group)
      } else {
        enough <- nrow(design_info$data) >= cfg$pseudobulk$min_samples_continuous
        group_counts <- NULL
      }
      fit <- if (enough) fit_linear(
        design_info$data$mean_score,
        design_info$design,
        design_info$coefficient
      ) else NULL
      if (is.null(fit)) {
        NULL
      } else {
        list(
          estimate = fit$estimate,
          standard_error = fit$se,
          conf_low = fit$low,
          conf_high = fit$high,
          statistic = fit$statistic,
          degrees_freedom = fit$df,
          p_value = fit$p,
          n_samples = nrow(design_info$data),
          n_reference = if (is.null(group_counts)) NA_integer_ else as.integer(group_counts[[reference_group]]),
          n_test = if (is.null(group_counts)) NA_integer_ else as.integer(group_counts[[test_group]]),
          covariates_used = paste(design_info$covariates_used, collapse = ";"),
          design_formula = design_info$formula
        )
      }
    }
  }, by = .(cell_type_analysis, program)]
  if (nrow(results) == 0L) return(data.table::data.table())
  results[, comparison := comparison]
  results[, effect := effect_column]
  results[, reference_group := reference_group %||% NA_character_]
  results[, test_group := test_group %||% NA_character_]
  results[, fdr := stats::p.adjust(p_value, method = "BH")]
  results[]
}

score_tests <- data.table::rbindlist(list(
  run_program_model(
    donor_scores, "Child_SLE_vs_HD", "condition", "binary",
    cfg$pseudobulk$covariates_child, "HD", "SLE", "Child"
  ),
  run_program_model(
    donor_scores, "Adult_SLE_vs_HD", "condition", "binary",
    cfg$pseudobulk$covariates_adult, "HD", "SLE", "Adult"
  ),
  run_program_model(
    donor_scores, "Child_SLE_High_vs_Low_Activity", "activity_binary", "binary",
    cfg$pseudobulk$covariates_activity, "Low", "High", "Child", "SLE"
  ),
  run_program_model(
    donor_scores, "Child_SLE_SLEDAI_continuous", "SLEDAI", "continuous",
    cfg$pseudobulk$covariates_activity, cohort = "Child", require_condition = "SLE"
  )
), use.names = TRUE, fill = TRUE)
if (nrow(score_tests) == 0L) {
  score_tests <- data.table::data.table(
    cell_type_analysis = character(), program = character(),
    estimate = numeric(), standard_error = numeric(), conf_low = numeric(),
    conf_high = numeric(), statistic = numeric(), degrees_freedom = numeric(),
    p_value = numeric(), n_samples = integer(), n_reference = integer(),
    n_test = integer(), covariates_used = character(), design_formula = character(),
    comparison = character(), effect = character(), reference_group = character(),
    test_group = character(), fdr = numeric()
  )
}
safe_write_tsv(score_tests, file.path(output_dir, "donor_level_program_models.tsv"))

child_heatmap_data <- score_tests[comparison == "Child_SLE_vs_HD" & !is.na(estimate)]
if (nrow(child_heatmap_data) > 0L) {
  child_heatmap <- ggplot2::ggplot(
    child_heatmap_data,
    ggplot2::aes(x = program, y = cell_type_analysis, fill = estimate)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.ticks = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Child cohort: adjusted donor-level SLE program effects",
      x = NULL, y = NULL, fill = "Adjusted effect"
    )
  ggplot2::ggsave(
    file.path(output_dir, "child_sle_vs_hd_program_score_heatmap.pdf"),
    child_heatmap, width = 11, height = 7
  )
}
if (!file.exists(file.path(output_dir, "child_sle_vs_hd_program_score_heatmap.pdf"))) {
  placeholder <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No estimable program-score effects")
  ggplot2::ggsave(
    file.path(output_dir, "child_sle_vs_hd_program_score_heatmap.pdf"),
    placeholder,
    width = 8,
    height = 4
  )
}

atomic_save_rds(object, output_file, compress = TRUE)
write_session_info(file.path(cfg$paths$logs, "06_program_scores_sessionInfo.txt"))
critical_outputs <- c(
  output_file,
  score_checkpoint,
  file.path(output_dir, "program_gene_matching.tsv"),
  file.path(output_dir, "donor_celltype_program_scores.tsv"),
  file.path(output_dir, "donor_level_program_models.tsv"),
  cell_exclusion_file,
  cell_exclusion_summary_file,
  file.path(output_dir, "child_sle_vs_hd_program_score_heatmap.pdf")
)
complete_step_and_cleanup(
  output_dir,
  step_signature,
  critical_outputs,
  cleanup_paths = c(
    file.path(output_dir, "program_column_map.tsv"),
    file.path(checkpoint_dir, "01_ucell_score_matrix.rds"),
    file.path(checkpoint_dir, "01_ucell_score_matrix.rds.signature.rds")
  )
)
log_message("SLE program scoring completed: ", output_file)
