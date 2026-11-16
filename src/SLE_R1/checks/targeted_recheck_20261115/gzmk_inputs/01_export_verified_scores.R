if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub('^--file=', '', grep('^--file=', args, value = TRUE)[1]), winslash = '/')
OUT <- dirname(script)
P <- 'F:/program/Hei/SLE/GSE121239_pipeline'
P135 <- 'F:/program/Hei/SLE/GSE135779_pipeline'
wt <- function(x, name) fwrite(x, file.path(OUT, name), sep = '\t', na = 'NA', quote = 'auto')
object <- readRDS(file.path(P, 'data/processed/01_gene_expression.rds'))
expression <- object$expression
meta <- as.data.table(object$metadata)
stopifnot(is.matrix(expression), identical(dim(expression), c(21914L, 312L)),
          !anyDuplicated(rownames(expression)), !anyDuplicated(colnames(expression)),
          !anyDuplicated(meta$sample_id), all(is.finite(expression)),
          setequal(meta$sample_id, colnames(expression)))
meta <- meta[match(colnames(expression), sample_id)]
modules <- fread('F:/program/Hei/SLE_v1/pipeline/3.5/data/processed/cross_cohort_module_membership.tsv')
labels <- c(cyan = 'CYAN_NEUTROPHIL_IMMATURE', greenyellow = 'GREENYELLOW_B_CELL_COMPLEMENT', magenta = 'MAGENTA_IFN_ANTIVIRAL')
primary <- unique(modules[module %in% names(labels), .(program = unname(labels[module]), gene = toupper(gene), module)])
primary[, `:=`(analysis_role = 'primary', source = 'GSE65391_frozen_WGCNA_module')]
secondary <- rbindlist(list(
  data.table(program = 'CD8_GZMK_STATE', gene = c('CD8A','CD8B','GZMK','CCL5','TRAC','CD3D','CD3E','LAG3')),
  data.table(program = 'CYTOTOXICITY_SHARED', gene = c('NKG7','GNLY','PRF1','GZMB','GZMH','CCL5','CTSW','KLRD1','KLRF1','FGFBP2'))))
secondary[, `:=`(module = NA_character_, analysis_role = 'secondary', source = 'GSE65391_predefined_cell_state')]
membership <- unique(rbindlist(list(primary, secondary), fill = TRUE))
membership[, matched := gene %in% rownames(expression)]
programs <- c(unname(labels), 'CYTOTOXICITY_SHARED', 'CD8_GZMK_STATE')
changed_programs <- c('CD8_GZMK_STATE', 'CYTOTOXICITY_SHARED')
ranks <- matrixStats::colRanks(expression, ties.method = 'average', preserveShape = TRUE, useNames = TRUE)
ranks <- ranks / (nrow(expression) + 1) - 0.5
dimnames(ranks) <- dimnames(expression)
scores <- t(vapply(programs, function(program_name) {
  genes <- membership[program == program_name & matched == TRUE, gene]
  colMeans(ranks[genes, , drop = FALSE])
}, numeric(ncol(expression))))
dimnames(scores) <- list(programs, colnames(expression))
no_ccl5 <- t(vapply(changed_programs, function(program_name) {
  genes <- membership[program == program_name & matched == TRUE & gene != 'CCL5', gene]
  colMeans(ranks[genes, , drop = FALSE])
}, numeric(ncol(expression))))
dimnames(no_ccl5) <- list(paste0(changed_programs, '_noCCL5'), colnames(expression))
coverage <- membership[, .(requested_genes = .N, matched_genes = sum(matched)), by = program]
coverage[, noCCL5_matched_genes := matched_genes - as.integer(program %in% changed_programs)]
stopifnot(coverage[program == 'CD8_GZMK_STATE', matched_genes] == 8,
          coverage[program == 'CYTOTOXICITY_SHARED', matched_genes] == 10,
          all(membership[program %in% changed_programs & gene == 'CCL5', matched]))
membership[, retained_in_noCCL5 := matched & !(program %in% changed_programs & gene == 'CCL5')]
membership[, original_weight := as.numeric(matched) / sum(matched), by = program]
membership[, noCCL5_weight := as.numeric(retained_in_noCCL5) / sum(retained_in_noCCL5), by = program]
wt(membership, 'program_membership.tsv')
wt(coverage, 'program_coverage.tsv')
wt(data.table(gene = rownames(expression), background_size = nrow(expression), rank_denominator = nrow(expression) + 1L), 'ranking_background_genes.tsv')
wt(as.data.table(object$selected_probes), 'original_selected_probes.tsv')
sle_meta <- meta[disease_state == 'Systemic Lupus Erythematosus']
stopifnot(nrow(sle_meta) == 292L, uniqueN(sle_meta$patient_id) == 65L,
          !anyNA(sle_meta$patient_id), !anyDuplicated(sle_meta$sample_id))
wt(sle_meta, 'sle_metadata_292_visits.tsv')
wt(sle_meta[, .(n_visits = .N, first_visit_order = min(visit_order), last_visit_order = max(visit_order)), by = patient_id], 'sle_patient_visit_counts.tsv')
wide <- as.data.table(t(rbind(scores, no_ccl5)[, sle_meta$sample_id, drop = FALSE]), keep.rownames = 'sample_id')
wide <- merge(sle_meta, wide, by = 'sample_id', sort = FALSE)
setorder(wide, patient_id, visit_order, sample_id)
wide[, CCL5_centered_fractional_rank := as.numeric(ranks['CCL5', sample_id])]
stopifnot(max(abs(wide$CD8_GZMK_STATE_noCCL5 - (8 * wide$CD8_GZMK_STATE - wide$CCL5_centered_fractional_rank) / 7)) < 1e-12,
          max(abs(wide$CYTOTOXICITY_SHARED_noCCL5 - (10 * wide$CYTOTOXICITY_SHARED - wide$CCL5_centered_fractional_rank) / 9)) < 1e-12)
wt(wide, 'sle_visit_scores.tsv')
s135 <- fread(file.path(P135, 'results/06_program_scores/donor_celltype_program_scores.tsv'))
child_sle <- s135[cohort == 'Child' & condition == 'SLE']
stopifnot(uniqueN(child_sle$subject_id) == 33L)
wt(child_sle, 'gse135779_child_sle_existing_scores.tsv')
wt(fread(file.path(P135, 'results/06_program_scores/program_gene_matching.tsv')), 'gse135779_program_gene_matching.tsv')
wt(fread(file.path(P135, 'config/sle_programs.tsv')), 'gse135779_requested_programs.tsv')
