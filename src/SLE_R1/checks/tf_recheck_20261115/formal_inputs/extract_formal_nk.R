if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")


suppressPackageStartupMessages({library(Matrix); library(data.table); library(edgeR)})
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub('^--file=', '', grep('^--file=', args, value = TRUE)[1]), winslash = '/')
out <- dirname(script)
cache <- 'F:/program/Hei/SLE/GSE135779_pipeline/data/processed/checkpoints/05_pseudobulk/pseudobulk_aggregations.rds'
sample_path <- 'F:/program/Hei/SLE/GSE135779_pipeline/data/metadata/sample_metadata.tsv'
members_path <- 'F:/program/Hei/SLE_v1/pipeline/3.5/data/processed/cross_cohort_module_membership.tsv'
supplement_path <- 'F:/program/Hei/SLE/metadata/frozen/2026-07-25_v1.0/GSE135779_published_supplement_extracted.tsv'
x <- readRDS(cache)$analysis
stopifnot(identical(colnames(x$counts), x$metadata$pseudobulk_id))
keep <- x$metadata$cohort == 'Child' & x$metadata$cell_label == 'NK'
md <- as.data.table(x$metadata[keep, , drop = FALSE])
counts <- as.matrix(x$counts[, keep, drop = FALSE])
stopifnot(nrow(md) == 44L, sum(md$condition == 'SLE') == 33L,
          sum(md$condition == 'HD') == 11L, sum(md$n_cells) == 15052L,
          !anyDuplicated(rownames(counts)), !anyDuplicated(md$subject_id),
          all(is.finite(counts)), all(counts >= 0), all(counts == round(counts)))

samples <- fread(sample_path, na.strings = c('NA', '', 'not_available', 'not_applicable'))
samples <- samples[match(md$sample_id, gsm)]
stopifnot(identical(samples$gsm, md$sample_id))
for (nm in c('SLEDAI', 'age')) {
  expected <- suppressWarnings(as.numeric(samples[[nm]]))
  observed <- as.numeric(md[[nm]])
  stopifnot(identical(is.na(expected), is.na(observed)), all(expected[!is.na(expected)] == observed[!is.na(observed)]))
}
md[, library_id := samples$library_id]
md[, raw_library_umi := colSums(counts)]
md[, mean_umi_per_cell := raw_library_umi / n_cells]
md[, log_mean_umi_per_cell := log(mean_umi_per_cell)]
md[, donor := subject_id]
colnames(counts) <- md$subject_id



clinical <- fread(supplement_path)[table == 'clinical']
clinical <- clinical[match(md$donor, Names)]
stopifnot(identical(clinical$Names, md$donor))
numeric_sledai <- suppressWarnings(as.numeric(clinical$SLEDAI))
sle <- md$condition == 'SLE'
stopifnot(identical(is.na(numeric_sledai[sle]), is.na(md$SLEDAI[sle])),
          all(numeric_sledai[sle & !is.na(numeric_sledai)] == md$SLEDAI[sle & !is.na(numeric_sledai)]),
          all(as.numeric(clinical$Age) == md$age),
          all(clinical$Batch == md$batch),
          all(ifelse(clinical$Gender == 'F', 'Female', 'Male') == md$sex))
fwrite(data.table(donor = md$donor, condition = md$condition,
                 clinical_sledai_literal = clinical$SLEDAI,
                 numeric_sledai = md$SLEDAI,
                 published_activity_group = md$published_activity_group,
                 analysis_activity_group = md$analysis_activity_group),
       file.path(out, 'formal_nk_clinical_provenance.tsv'), sep = '\t', quote = FALSE, na = 'NA')



nonzero <- rowSums(counts) > 0
y_for_norm <- DGEList(counts = counts[nonzero, , drop = FALSE])
y_for_norm <- calcNormFactors(y_for_norm, method = 'TMM')
y <- DGEList(counts = counts)
y$samples$norm.factors <- y_for_norm$samples$norm.factors
stopifnot(isTRUE(all.equal(y$samples$lib.size, md$raw_library_umi)))
lcpm <- cpm(y, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
raw_cpm <- cpm(counts, log = FALSE)
md[, tmm_norm_factor := y$samples$norm.factors]
md[, effective_library_umi := raw_library_umi * tmm_norm_factor]

write_mat <- function(mat, name) {
  d <- as.data.table(mat)
  d[, gene := rownames(mat)]
  setcolorder(d, c('gene', colnames(mat)))
  fwrite(d, file.path(out, name), sep = '\t', quote = FALSE, na = 'NA')
}
write_mat(counts, 'formal_nk_counts.tsv')
write_mat(lcpm, 'formal_nk_tmm_logcpm.tsv')
fwrite(md, file.path(out, 'formal_nk_metadata.tsv'), sep = '\t', quote = FALSE, na = 'NA')
genes <- data.table(gene = rownames(counts), total_counts = rowSums(counts),
                    donors_detected = rowSums(counts > 0),
                    donors_raw_cpm_ge1 = rowSums(raw_cpm >= 1),
                    eligible_cpm1_in5 = rowSums(raw_cpm >= 1) >= 5)
fwrite(genes, file.path(out, 'formal_nk_gene_detection.tsv'), sep = '\t', quote = FALSE, na = 'NA')
members <- fread(members_path)
ifn <- members[module == 'magenta', gene]
stopifnot(length(ifn) == 157L, !anyDuplicated(ifn))
fwrite(data.table(gene = ifn, present_in_formal_counts = ifn %in% rownames(counts)),
       file.path(out, 'formal_frozen_ifn_members.tsv'), sep = '\t', quote = FALSE)

input_paths <- c(cache, sample_path, members_path, supplement_path, script)
manifest <- data.table(path = input_paths, sha256 = vapply(input_paths, function(p) digest::digest(file = p, algo = 'sha256'), character(1)))
fwrite(manifest, file.path(out, 'input_sha256.tsv'), sep = '\t', quote = FALSE)
settings <- data.table(parameter = c('selection', 'source_level', 'donors', 'SLE', 'HD', 'NK_cells', 'genes',
                               'tmm_method', 'tmm_gene_set', 'logCPM_prior_count', 'logCPM_base',
                               'mean_umi_definition', 'SLEDAI_missing_SLE', 'SLEDAI_missing_HD', 'edgeR_version'),
                       value = c('cohort Child AND cell_label NK', 'analysis', nrow(md), sum(md$condition == 'SLE'),
                                 sum(md$condition == 'HD'), sum(md$n_cells), nrow(counts), 'TMM',
                                 paste0(sum(nonzero), ' genes with nonzero total counts'), 2, 2,
                                 'raw column sum / n_cells', sum(is.na(md$SLEDAI) & md$condition == 'SLE'),
                                 sum(is.na(md$SLEDAI) & md$condition == 'HD'), as.character(packageVersion('edgeR'))))
fwrite(settings, file.path(out, 'normalization_settings.tsv'), sep = '\t', quote = FALSE)
capture.output(sessionInfo(), file = file.path(out, 'R_sessionInfo.txt'))
cat('Exported', nrow(counts), 'genes x', ncol(counts), 'donors;', sum(md$n_cells), 'formal QC NK cells.\n')
cat('SLE numeric SLEDAI:', sum(md$condition == 'SLE' & !is.na(md$SLEDAI)), 'of', sum(md$condition == 'SLE'), '\n')
cat('Frozen IFN matching:', sum(ifn %in% rownames(counts)), '/', length(ifn), '\n')
