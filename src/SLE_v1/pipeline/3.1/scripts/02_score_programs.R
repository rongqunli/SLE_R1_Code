if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(c("data.table", "matrixStats", "ggplot2"), "frozen program scoring")
source(file.path(PIPELINE_ROOT, "scripts", "_program_helpers.R"), local = TRUE)

script_path <- script_file()
expression_file <- file.path(cfg$paths$processed, "01_gene_expression.rds")
program_file <- cfg$files$programs
output_dir <- file.path(cfg$paths$results, "02_program_scores")
object_file <- file.path(cfg$paths$processed, "02_program_scores.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- c(
  object_file,
  file.path(output_dir, c(
    "program_gene_matching.tsv", "sample_program_scores.tsv",
    "program_score_summary.tsv", "program_score_distributions.pdf"
  ))
)
signature <- build_step_signature(
  input_files = c(expression_file, program_file),
  code_files = c(script_path, file.path(PIPELINE_ROOT, "scripts", "_common.R"),
                 file.path(PIPELINE_ROOT, "scripts", "_program_helpers.R"),
                 file.path(PIPELINE_ROOT, "config", "config.R"))
)
if (!has_flag("--force") && step_is_complete(output_dir, signature, outputs)) {
  log_message("Step 02 frozen program scoring is already complete.")
  quit(save = "no", status = 0L)
}
if (!file.exists(expression_file)) stop("Missing step 01 expression object.")

expression_object <- readRDS(expression_file)
expression <- expression_object$expression
metadata <- data.table::as.data.table(expression_object$metadata)
selected_probes <- data.table::as.data.table(expression_object$selected_probes)
programs <- data.table::fread(program_file, sep = "\t")
if (anyDuplicated(programs$program)) stop("Duplicate program names in config.")

matching <- data.table::rbindlist(lapply(seq_len(nrow(programs)), function(i) {
  configured <- unique(toupper(trimws(strsplit(programs$genes[[i]], ";", fixed = TRUE)[[1]])))
  data.table::data.table(
    program = programs$program[[i]],
    category = programs$category[[i]],
    primary = programs$primary[[i]],
    gene = configured,
    matched = configured %in% rownames(expression)
  )
}))
matching <- merge(
  matching,
  selected_probes[, .(gene = gene_symbol, probe_id, probe_variance)],
  by = "gene", all.x = TRUE, sort = FALSE
)
data.table::setcolorder(matching, c(
  "program", "category", "primary", "gene", "matched", "probe_id", "probe_variance"
))
coverage <- matching[, .(matched_genes = sum(matched), configured_genes = .N), by = program]
if (any(coverage$matched_genes < cfg$expression$min_program_genes)) {
  failed <- coverage[matched_genes < cfg$expression$min_program_genes, program]
  stop("Insufficient matched genes for program(s): ", paste(failed, collapse = ", "))
}
safe_write_tsv(matching, file.path(output_dir, "program_gene_matching.tsv"))

log_message("Ranking genes within each of ", ncol(expression), " included arrays.")
ranks <- matrixStats::colRanks(
  expression,
  ties.method = "average",
  preserveShape = TRUE,
  useNames = TRUE
)
ranks <- ranks / (nrow(expression) + 1) - 0.5
rownames(ranks) <- rownames(expression)
colnames(ranks) <- colnames(expression)

score_matrix <- vapply(programs$program, function(program) {
  genes <- matched_program_genes(matching, program)
  colMeans(ranks[genes, , drop = FALSE])
}, numeric(ncol(expression)))
score_matrix <- t(score_matrix)
rownames(score_matrix) <- programs$program
colnames(score_matrix) <- colnames(expression)
rm(ranks); gc(FALSE)

score_long <- data.table::as.data.table(as.table(score_matrix))
data.table::setnames(score_long, c("program", "sample_id", "score"))
score_long[, `:=`(program = as.character(program), sample_id = as.character(sample_id))]
score_long <- merge(score_long, metadata, by = "sample_id", all.x = TRUE, sort = FALSE)
score_long <- score_long[match(
  paste(rep(programs$program, each = ncol(score_matrix)), rep(colnames(score_matrix), times = nrow(score_matrix))),
  paste(program, sample_id)
)]
if (anyNA(score_long$patient_id)) stop("Program scores failed to join frozen metadata.")

score_summary <- score_long[, .(
  n_arrays = .N,
  mean_score = mean(score),
  sd_score = stats::sd(score),
  median_score = stats::median(score),
  min_score = min(score),
  max_score = max(score)
), by = .(program, disease_state)]

safe_write_tsv(score_long, file.path(output_dir, "sample_program_scores.tsv"))
safe_write_tsv(score_summary, file.path(output_dir, "program_score_summary.tsv"))

pdf(file.path(output_dir, "program_score_distributions.pdf"), width = 11, height = 8.5, onefile = TRUE)
print(ggplot2::ggplot(score_long, ggplot2::aes(score, fill = disease_state)) +
  ggplot2::geom_density(alpha = 0.35) +
  ggplot2::facet_wrap(~ program, scales = "free", ncol = 3) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Frozen rank-based program scores", fill = "State"))
print(ggplot2::ggplot(score_long, ggplot2::aes(disease_state, score, fill = disease_state)) +
  ggplot2::geom_boxplot(outlier.size = 0.4) +
  ggplot2::facet_wrap(~ program, scales = "free", ncol = 3) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
  ggplot2::labs(title = "Program score distributions by disease state", x = NULL, fill = "State"))
grDevices::dev.off()

score_object <- list(
  score_matrix = score_matrix,
  metadata = as.data.frame(metadata),
  matching = as.data.frame(matching),
  method = "within_sample_fractional_rank_mean_centered_at_0.5"
)
atomic_save_rds(score_object, object_file, compress = TRUE)
write_session_info(file.path(cfg$paths$logs, "02_program_scores_sessionInfo.txt"))
mark_step_complete(output_dir, signature, outputs)
log_message("Program scoring complete: ", nrow(programs), " programs x ",
            ncol(score_matrix), " included arrays.")
