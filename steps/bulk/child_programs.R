library(data.table)
root <- Sys.getenv('SLE_RUN_ROOT')
stopifnot(nzchar(root))
membership <- fread(file.path(root, 'SLE_v1/results/3.1/06_wgcna/module_membership.tsv'))
stopifnot(!anyDuplicated(membership$gene))
program_map <- data.table(program = c('MAGENTA_IFN_ANTIVIRAL', 'GREENYELLOW_B_CELL_COMPLEMENT', 'CYAN_NEUTROPHIL_IMMATURE'), module = c('magenta', 'greenyellow', 'cyan'))
expanded <- merge(program_map, membership[module %in% program_map$module, .(gene, module, kME, hub_rank)], by = 'module')
setorder(expanded, program, -kME, gene)
for (relative in c('SLE/GSE200306_pipeline/results/07_cross_dataset_bridge/expanded_program_membership.tsv', 'SLE_v1/pipeline/3.2/data/reference/expanded_program_membership.tsv')) {
    target <- file.path(root, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    fwrite(expanded, target, sep = '\t')
}
for (pair in list(c('01_gene_expression.rds', 'GSE65391_gene_expression.rds'), c('04_wgcna.rds', 'GSE65391_wgcna.rds'))) {
    source <- file.path(root, 'SLE_v1/pipeline/3.1/data/processed', pair[[1]])
    target <- file.path(root, 'SLE_v1/pipeline/3.2/data/reference', pair[[2]])
    stopifnot(file.copy(source, target, overwrite = TRUE))
}
