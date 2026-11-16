if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
options(stringsAsFactors=FALSE,scipen=99)
suppressPackageStartupMessages({library(data.table);library(matrixStats);library(limma);library(WGCNA);library(fgsea)})
project_root <- "F:/program/Hei/SLE_v1"
table_dir <- file.path(project_root,"Figures/tables")
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
target_modules <- c("magenta","greenyellow","cyan")
seed <- 20260726L
read_project_tsv <- function(path) fread(file.path(project_root,path))
center_expression_by_patient <- function(expression, patient_id) {
  centered <- expression
  patient_columns <- split(seq_along(patient_id), patient_id)
  for (columns in patient_columns) {
    centered[, columns] <- sweep(
      expression[, columns, drop = FALSE], 1L,
      rowMeans(expression[, columns, drop = FALSE]), FUN = "-"
    )
  }
  centered
}

message("Reconstructing the exact WGCNA input matrix.")
expression_object <- readRDS(file.path(
  project_root, "pipeline/3.1/data/processed/01_gene_expression.rds"
))
expression <- expression_object$expression
clinical <- read_project_tsv("results/3.1/03_longitudinal/longitudinal_analysis_data.tsv")
repeated_data <- copy(clinical[patient_visit_count >= 2L])
columns <- match(repeated_data$sample_id, colnames(expression))
stopifnot(nrow(repeated_data) == 897L, uniqueN(repeated_data$patient_id) == 131L,
          !anyNA(columns))
repeated_expression <- expression[, columns, drop = FALSE]
repeated_data[, `:=`(
  patient_id = factor(patient_id),
  batch = factor(as.character(batch))
)]
centered <- center_expression_by_patient(repeated_expression, repeated_data$patient_id)
batch_residualized <- limma::removeBatchEffect(centered, batch = repeated_data$batch)
batch_residualized <- center_expression_by_patient(
  batch_residualized, repeated_data$patient_id
)
selection_variance <- matrixStats::rowVars(batch_residualized)
selected_indices <- head(order(selection_variance, decreasing = TRUE), 5000L)
selected_expression <- batch_residualized[selected_indices, , drop = FALSE]
dat_expr <- t(selected_expression)
rownames(dat_expr) <- repeated_data$sample_id
stopifnot(WGCNA::goodSamplesGenes(dat_expr, verbose = 0)$allOK)

assignments <- read_project_tsv("results/3.1/06_wgcna/module_assignments.tsv")
soft_diagnostics <- read_project_tsv("results/3.1/06_wgcna/soft_threshold_diagnostics.tsv")
associations <- read_project_tsv("results/3.1/06_wgcna/module_trait_associations.tsv")
stability <- read_project_tsv("results/3.1/06_wgcna/module_stability.tsv")
membership <- read_project_tsv("results/3.1/06_wgcna/module_membership.tsv")
hubs <- read_project_tsv("results/3.1/06_wgcna/hub_genes.tsv")
chosen_power <- soft_diagnostics[selected == TRUE, power][[1L]]

go_pathways <- fgsea::gmtPathways(file.path(
  project_root, "pipeline/3.2/config/gene_sets/c5.go.symbols.gmt"
))
kegg_pathways <- fgsea::gmtPathways(file.path(
  project_root, "pipeline/3.2/config/gene_sets/c2.cp.kegg.symbols.gmt"
))
pathway_collections <- list(GO = go_pathways, KEGG = kegg_pathways)
universe <- unique(assignments$gene)

ontology_from_pathway <- function(pathway) {
  data.table::fcase(
    grepl("^GOBP_", pathway), "GO BP",
    grepl("^GOCC_", pathway), "GO CC",
    grepl("^GOMF_", pathway), "GO MF",
    grepl("^KEGG_", pathway), "KEGG",
    default = "Other"
  )
}

run_ora <- function(module_name, collection_name, pathways) {
  module_genes <- assignments[module == module_name, gene]
  result <- rbindlist(lapply(names(pathways), function(pathway_name) {
    pathway_genes <- intersect(unique(pathways[[pathway_name]]), universe)
    overlap <- intersect(module_genes, pathway_genes)
    K <- length(module_genes)
    M <- length(pathway_genes)
    N <- length(universe)
    k <- length(overlap)
    if (M < 5L || M > 1000L || k < 1L) return(NULL)
    background_other <- N - M - K + k
    odds_haldane <- ((k + 0.5) * (background_other + 0.5)) /
      ((K - k + 0.5) * (M - k + 0.5))
    data.table(
      module = module_name, collection = collection_name,
      pathway = pathway_name, module_size = K, pathway_size_in_universe = M,
      overlap = k, gene_ratio = k / K, background_ratio = M / N,
      fold_enrichment = (k / K) / (M / N),
      odds_ratio_haldane = odds_haldane,
      p_value = phyper(k - 1L, M, N - M, K, lower.tail = FALSE),
      overlap_genes = paste(overlap, collapse = ";")
    )
  }), fill = TRUE)
  result[, fdr := p.adjust(p_value, method = "BH")]
  setorder(result, fdr, -fold_enrichment, -overlap)
  result
}

ora_file <- file.path(table_dir, "selected_modules_GO_KEGG_overrepresentation.tsv")
if (file.exists(ora_file)) {
  ora <- fread(ora_file, sep = "\t")
} else {
  ora <- rbindlist(lapply(target_modules, function(module_name) {
    rbindlist(lapply(names(pathway_collections), function(collection_name) {
      run_ora(module_name, collection_name, pathway_collections[[collection_name]])
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
  fwrite(ora, ora_file, sep = "\t")
}
ora[, ontology := ontology_from_pathway(pathway)]
setcolorder(ora, c("module", "collection", "ontology",
                   setdiff(names(ora), c("module", "collection", "ontology"))))
fwrite(ora, ora_file, sep = "\t")

message("Preparing module-membership-ranked GO and KEGG GSEA.")

run_module_gsea <- function(module_name, collection_name, pathways) {
  eigengene_column <- paste0("ME", module_name)
  ranks <- stats::cor(
    dat_expr, eigengene_table[[eigengene_column]][eigengene_rows],
    use = "pairwise.complete.obs"
  )[, 1L]
  names(ranks) <- colnames(dat_expr)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  set.seed(seed + match(module_name, target_modules))
  ans <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = pathways, stats = ranks, minSize = 10L,
    maxSize = 500L, eps = 1e-20
  ))
  ans[, `:=`(module = module_name, collection = collection_name)]
  ans[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1L))]
  ans[]
}

gsea_file <- file.path(table_dir, "selected_modules_kME_ranked_GO_KEGG_GSEA.tsv")
if (file.exists(gsea_file)) {
  gsea <- fread(gsea_file, sep = "\t")
} else {
  wgcna_object <- readRDS(file.path(
    project_root, "pipeline/3.1/data/processed/04_wgcna.rds"
  ))
  eigengene_table <- as.data.table(wgcna_object$module_eigengenes)
  eigengene_rows <- match(rownames(dat_expr), eigengene_table$sample_id)
  stopifnot(!anyNA(eigengene_rows))
  gsea <- rbindlist(lapply(target_modules, function(module_name) {
    rbindlist(lapply(names(pathway_collections), function(collection_name) {
      run_module_gsea(module_name, collection_name, pathway_collections[[collection_name]])
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
  setcolorder(gsea, c("module", "collection", setdiff(names(gsea), c("module", "collection"))))
  setorder(gsea, module, collection, padj, -NES)
  fwrite(gsea, gsea_file, sep = "\t")
}
gsea[, ontology := ontology_from_pathway(pathway)]
setcolorder(gsea, c("module", "collection", "ontology",
                    setdiff(names(gsea), c("module", "collection", "ontology"))))
fwrite(gsea, gsea_file, sep = "\t")

clean_term <- function(x) {
  x <- sub("^(GOBP|GOCC|GOMF|KEGG)_", "", x)
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("B Cell", "B-cell", x, fixed = TRUE)
  x <- gsub("b Cell", "B-cell", x, fixed = TRUE)
  x <- gsub("T Cell", "T-cell", x, fixed = TRUE)
  x <- gsub("t Cell", "T-cell", x, fixed = TRUE)
  x <- gsub("Type i ", "Type I ", x, fixed = TRUE)
  x <- gsub("Interferon Alpha", "Interferon-alpha", x, fixed = TRUE)
  x <- gsub("Interferon Beta", "Interferon-beta", x, fixed = TRUE)
  x <- gsub("Interferon Mediated", "Interferon-mediated", x, fixed = TRUE)
  x <- gsub("Rig i Like Receptor", "RIG-I-like receptor", x, fixed = TRUE)
  x <- gsub("Rna", "RNA", x, fixed = TRUE)
  x <- gsub("Dna", "DNA", x, fixed = TRUE)
  x
}

representative_pathways <- list(
  magenta = c(
    "GOBP_RESPONSE_TO_VIRUS",
    "GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY",
    "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
    "GOBP_NEGATIVE_REGULATION_OF_VIRAL_GENOME_REPLICATION",
    "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
    "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY"
  ),
  greenyellow = c(
    "GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_B_CELL_ACTIVATION",
    "GOBP_B_CELL_DIFFERENTIATION",
    "GOBP_HUMORAL_IMMUNE_RESPONSE",
    "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
    "KEGG_HEMATOPOIETIC_CELL_LINEAGE"
  ),
  cyan = c(
    "GOCC_SECRETORY_GRANULE",
    "GOCC_SPECIFIC_GRANULE",
    "GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE",
    "GOBP_DEFENSE_RESPONSE_TO_BACTERIUM",
    "GOBP_NEUTROPHIL_MEDIATED_CYTOTOXICITY",
    "GOBP_NEUTROPHIL_MEDIATED_IMMUNITY"
  )
)
annotation_terms <- rbindlist(lapply(target_modules, function(module_name) {
  selected <- ora[
    module == module_name & pathway %chin% representative_pathways[[module_name]] &
      fdr < 0.05 & overlap >= 3L
  ]
  selected[, preferred_order := match(pathway, representative_pathways[[module_name]])]
  setorder(selected, preferred_order)
  selected
}), use.names = TRUE, fill = TRUE)
annotation_terms[, display_term := paste0(clean_term(pathway), " (", ontology, ")")]
annotation_terms[, neglog10_fdr := -log10(pmax(fdr, 1e-300))]
annotation_terms[, fdr_label := paste0(
  "FDR=", formatC(fdr, format = "e", digits = 1)
)]
annotation_terms[, selection_rule := paste0(
  "FDR<0.05; overlap>=3; representative nonredundant hypothesis-aligned terms; full results retained"
)]
fwrite(annotation_terms, file.path(table_dir, "selected_annotation_terms_for_Figure1.tsv"), sep = "\t")
canonical_dir <- file.path(project_root,"results/3.1/07_wgcna_functional_annotation")
dir.create(canonical_dir,recursive=TRUE,showWarnings=FALSE)
for (name in c("selected_modules_GO_KEGG_overrepresentation.tsv","selected_modules_kME_ranked_GO_KEGG_GSEA.tsv","selected_annotation_terms_for_Figure1.tsv")) {
  stopifnot(file.copy(file.path(table_dir,name),file.path(canonical_dir,name),overwrite=TRUE))
}

