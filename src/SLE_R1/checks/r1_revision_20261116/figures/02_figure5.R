if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")

options(stringsAsFactors=FALSE)
stage <- 'F:/program/Hei/SLE_R1/checks/r1_revision_20261116'
dir.create(file.path(stage,'work'),recursive=TRUE,showWarnings=FALSE)
options(stringsAsFactors = FALSE, scipen = 99)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrastr)
  library(ggalluvial)
  library(ggforce)
  library(ggrepel)
  library(viridis)
})

args <- c('F:/program/Hei/SLE_v1',file.path(stage,'figures/Figure5_panels'),file.path(stage,'figures/Figure5_records'))
project_root <- if (length(args) >= 1L) args[[1L]] else "F:/program/Hei/SLE_v1"
output_root <- if (length(args) >= 2L) args[[2L]] else file.path(project_root, "Figures", "Figure5")
record_root <- if (length(args) >= 3L) args[[3L]] else file.path(
  project_root, "results", "3.5", "07_figure5_group_context"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(record_root, recursive = TRUE, showWarnings = FALSE)

font_family <- "Times New Roman"
condition_colors <- c(HD = "#2C6E9B", SLE = "#C23B3B")
program_labels <- c(
  magenta = "IFN/antiviral",
  greenyellow = "B-cell/BCR",
  cyan = "Neutrophil granule/antimicrobial"
)

theme_article <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, colour = "#202124"),
      plot.title = element_text(face = "bold", size = rel(1.16), margin = margin(b = 7)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.text = element_text(colour = "#202124"),
      axis.title = element_text(face = "plain"),
      strip.background = element_rect(fill = "#F3F4F5", colour = "#D5D8DC", linewidth = 0.45),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(9, 11, 9, 9)
    )
}

theme_nebula <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, colour = "white"),
      panel.background = element_rect(fill = "black", colour = NA),
      panel.border = element_rect(fill = NA, colour = "white", linewidth = 0.35),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "black", colour = "black"),
      plot.title = element_text(face = "bold", colour = "white", size = rel(1.16), margin = margin(b = 2)),
      plot.subtitle = element_text(colour = "#CFCFCF", size = rel(0.88), margin = margin(b = 7)),
      axis.title = element_text(colour = "white"),
      axis.text = element_text(colour = "white"),
      axis.ticks = element_line(colour = "white", linewidth = 0.25),
      strip.background = element_rect(fill = "#272727", colour = "#777777", linewidth = 0.35),
      strip.text = element_text(face = "bold", colour = "white", lineheight = 0.95, margin = margin(6, 6, 6, 6)),
      legend.position = "none",
      panel.spacing = grid::unit(0.8, "lines"),
      plot.margin = margin(9, 9, 9, 9),
      plot.tag = element_text(family = font_family, face = "bold", colour = "white")
    )
}

theme_set(theme_article())
update_geom_defaults("text", list(family = font_family, colour = "#202124"))
update_geom_defaults("label", list(family = font_family, colour = "#202124"))

manifest <- list()
noop_panel <- function(...) invisible(NULL)

read_project_tsv <- function(path) {
  fread(file.path(project_root, path), sep = "\t", na.strings = c("", "NA"))
}
fmt_p <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x))

safe_bandwidth <- function(x, adjust = 1) {
  bandwidth <- suppressWarnings(MASS::bandwidth.nrd(x))
  data_range <- diff(range(x, finite = TRUE))
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    bandwidth <- max(data_range / 30, 1e-3)
  }
  bandwidth * adjust
}

density_grid_one <- function(data, x_limits, y_limits, n = 260, adjust = 1.20) {
  if (nrow(data) < 3L) {
    return(data.table(UMAP_1 = numeric(), UMAP_2 = numeric(), density = numeric()))
  }
  kde <- MASS::kde2d(
    x = data$UMAP_1,
    y = data$UMAP_2,
    h = c(
      safe_bandwidth(data$UMAP_1, adjust),
      safe_bandwidth(data$UMAP_2, adjust)
    ),
    n = n,
    lims = c(x_limits, y_limits)
  )
  grid <- as.data.table(expand.grid(UMAP_1 = kde$x, UMAP_2 = kde$y))
  grid[, density := as.vector(kde$z)]
  grid
}

full_umap_path <- "results/exploratory_program_landscape/05_single_cell_trajectory/checkpoints/full_atlas_umap_metadata.rds"
full_umap_file <- file.path(project_root, full_umap_path)
full_umap <- as.data.table(readRDS(full_umap_file))
required_columns <- c(
  "cell_id", "UMAP_1", "UMAP_2", "sample_id", "subject_id",
  "cohort", "condition", "cell_type_analysis"
)
if (!all(required_columns %chin% names(full_umap))) {
  stop("Frozen atlas metadata is missing: ", paste(setdiff(required_columns, names(full_umap)), collapse = ", "))
}
if (anyDuplicated(full_umap$cell_id)) stop("Frozen atlas contains duplicated cell_id values")
if (any(!is.finite(full_umap$UMAP_1)) || any(!is.finite(full_umap$UMAP_2))) {
  stop("Frozen atlas contains non-finite UMAP coordinates")
}

cell_type_order <- c(
  "B", "B_Naive", "Plasmablast", "Monocyte", "CD14_Monocyte", "CD16_Monocyte", "pDC",
  "CD4_T", "CD4_T_Naive", "CD8_T_Naive", "CD8_T_GZMK", "CD8_T_Cytotoxic",
  "NK", "Platelet", "Erythroid"
)
cell_type_labels <- c(
  B = "B", B_Naive = "Naive B", Plasmablast = "Plasmablast", Monocyte = "Monocyte",
  CD14_Monocyte = "CD14 monocyte", CD16_Monocyte = "CD16 monocyte", pDC = "pDC",
  CD4_T = "CD4 T", CD4_T_Naive = "Naive CD4 T", CD8_T_Naive = "Naive CD8 T",
  CD8_T_GZMK = "GZMK CD8 T", CD8_T_Cytotoxic = "Cytotoxic CD8 T", NK = "NK",
  Platelet = "Platelet", Erythroid = "Erythroid"
)
cell_type_short <- c(
  B = "B", B_Naive = "Naive B", Plasmablast = "PB", Monocyte = "Mono",
  CD14_Monocyte = "CD14", CD16_Monocyte = "CD16", pDC = "pDC",
  CD4_T = "CD4 T", CD4_T_Naive = "Naive CD4", CD8_T_Naive = "Naive CD8",
  CD8_T_GZMK = "GZMK CD8", CD8_T_Cytotoxic = "Cyto CD8", NK = "NK",
  Platelet = "Plt", Erythroid = "Ery"
)
cell_palette <- setNames(
  c(
    "#4E79A7", "#76B7B2", "#59A14F", "#FF9DA7", "#E15759", "#FF9D9A", "#B07AA1",
    "#F28E2B", "#EDC948", "#9C755F", "#BAB0AC", "#7F6AAD", "#2F4B7C",
    "#D37295", "#8CD17D"
  ),
  cell_type_order
)
missing_types <- setdiff(unique(full_umap$cell_type_analysis), cell_type_order)
if (length(missing_types)) stop("Unrecognized frozen cell type(s): ", paste(missing_types, collapse = ", "))

full_umap[, condition := factor(as.character(condition), levels = c("HD", "SLE"))]
full_umap[, cohort := factor(as.character(cohort), levels = c("Child", "Adult"))]
full_umap[, cell_type_label := factor(
  unname(cell_type_labels[cell_type_analysis]),
  levels = unname(cell_type_labels[cell_type_order])
)]
full_umap[, cell_type_short := factor(
  unname(cell_type_short[cell_type_analysis]),
  levels = unname(cell_type_short[cell_type_order])
)]

shared_x_limits <- range(full_umap$UMAP_1)
shared_y_limits <- range(full_umap$UMAP_2)
shared_x_scale <- scale_x_continuous(limits = shared_x_limits, expand = expansion(mult = 0.025))
shared_y_scale <- scale_y_continuous(limits = shared_y_limits, expand = expansion(mult = 0.025))

p5a <- ggplot(full_umap, aes(UMAP_1, UMAP_2, colour = cell_type_label)) +
  ggrastr::geom_point_rast(size = 0.24, alpha = 0.58, raster.dpi = 360) +
  scale_colour_manual(
    values = setNames(cell_palette, unname(cell_type_labels[names(cell_palette)])),
    name = NULL, na.value = "#BBBBBB"
  ) +
  shared_x_scale + shared_y_scale +
  labs(title = "Pediatric and adult PBMC atlas", x = "UMAP 1", y = "UMAP 2") +
  theme_article(9.5) +
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.key.width = grid::unit(0.9, "lines")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.7, alpha = 1)))
noop_panel(
  p5a, "A", "Figure5A_full_single_cell_atlas", "Pediatric and adult PBMC atlas",
  full_umap_path, 10.8, 7.2
)

group_counts <- full_umap[, .(cells = .N, donors = uniqueN(subject_id)), by = condition]
setorder(group_counts, condition)
group_counts[, condition_label_text := sprintf(
  "%s\n(n = %s cells; %d donors)", condition, comma(cells), donors
)]
group_counts[, condition_label := factor(condition_label_text, levels = condition_label_text)]
group_counts[, condition_label_text := NULL]
full_umap <- merge(full_umap, group_counts[, .(condition, condition_label)], by = "condition", all.x = TRUE, sort = FALSE)





global_box_area <- diff(shared_x_limits) * diff(shared_y_limits)
cell_type_counts <- full_umap[, .(n_total = .N), by = .(condition_label, cell_type_short)]
trimmed_hull_data <- full_umap[, {
  distance <- (UMAP_1 - median(UMAP_1))^2 + (UMAP_2 - median(UMAP_2))^2
  keep_n <- max(10L, ceiling(.N * 0.85))
  .SD[order(distance)[seq_len(keep_n)]]
}, by = .(condition_label, cell_type_short)]
hull_metrics <- trimmed_hull_data[, .(
  robust_box_area = diff(range(UMAP_1)) * diff(range(UMAP_2))
), by = .(condition_label, cell_type_short)]
hull_metrics <- merge(
  hull_metrics, cell_type_counts,
  by = c("condition_label", "cell_type_short"), all.x = TRUE, sort = FALSE
)
valid_hulls <- hull_metrics[n_total >= 1000 & robust_box_area / global_box_area <= 0.18,
  .(condition_label, cell_type_short)]
hull_data <- merge(
  trimmed_hull_data, valid_hulls,
  by = c("condition_label", "cell_type_short"), all = FALSE, sort = FALSE
)[, {
  if (.N > 1800L) .SD[unique(round(seq(1, .N, length.out = 1800L)))] else .SD
}, by = .(condition_label, cell_type_short)]
label_source <- merge(
  full_umap, cell_type_counts[n_total >= 100],
  by = c("condition_label", "cell_type_short"), all = FALSE, sort = FALSE
)
label_data <- label_source[, {
  distance <- (UMAP_1 - median(UMAP_1))^2 + (UMAP_2 - median(UMAP_2))^2
  keep_n <- max(1L, ceiling(.N * 0.75))
  core <- .SD[order(distance)[seq_len(keep_n)]]
  .(label_x = median(core$UMAP_1), label_y = median(core$UMAP_2))
}, by = .(condition_label, cell_type_short)]
rm(trimmed_hull_data, label_source); invisible(gc())

density_levels <- levels(full_umap$condition_label)
density_data <- rbindlist(lapply(density_levels, function(group_name) {
  group_data <- full_umap[as.character(condition_label) == group_name]
  grid <- density_grid_one(
    group_data,
    x_limits = shared_x_limits,
    y_limits = shared_y_limits,
    n = 260,
    adjust = 1.20
  )
  grid[, condition_label := factor(group_name, levels = density_levels)]
  grid
}))
density_limit <- max(density_data$density, na.rm = TRUE)
bright_core_labels <- c(
  "CD14", "NK", "Cyto CD8", "CD4 T", "Naive CD8", "Naive CD4", "Naive B"
)
label_data[, label_colour := fifelse(
  cell_type_short %chin% bright_core_labels, "black", "white"
)]
set.seed(20260810)
nebula_points <- full_umap[sample.int(nrow(full_umap))]

p5b <- ggplot() +
  geom_raster(
    data = density_data,
    aes(UMAP_1, UMAP_2, fill = density),
    interpolate = TRUE
  ) +
  ggrastr::geom_point_rast(
    data = nebula_points, aes(UMAP_1, UMAP_2),
    colour = "white", size = 0.012, alpha = 0.12, raster.dpi = 360
  ) +
  ggforce::geom_mark_hull(
    data = hull_data,
    aes(UMAP_1, UMAP_2, group = cell_type_short),
    linewidth = 0.42, concavity = 2, colour = scales::alpha("white", 0.82),
    expand = grid::unit(1.2, "mm"), linetype = "dashed",
    fill = NA, show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(label_x, label_y, label = cell_type_short, colour = label_colour),
    size = 2.55, fontface = "bold",
    box.padding = 0.34, point.padding = 0.10,
    segment.colour = scales::alpha("#B8B8B8", 0.62), segment.size = 0.28,
    min.segment.length = 0, max.overlaps = Inf, seed = 20260718,
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  viridis::scale_fill_viridis(
    option = "magma", direction = 1, limits = c(0, density_limit),
    oob = squish, guide = "none"
  ) +
  shared_x_scale + shared_y_scale +
  facet_wrap(~condition_label, nrow = 1) +
  labs(
    title = "HD and SLE on the frozen PBMC atlas",
    subtitle = "Shared density scale and frozen UMAP coordinates",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_nebula(9.5)
p5b <- p5b + theme(plot.tag = element_text(colour = "white"))
noop_panel(
  p5b, "B", "Figure5B_HD_SLE_nebula_atlas",
  "HD and SLE on the frozen PBMC atlas",
  c(full_umap_path, "reference: Diabetes2/16.T2D/snRNA/12_publication_ready_figures/06_myonuclei_Nebula_by_RiskGroup.pdf"),
  13.2, 5.4, background = "black"
)

donor_meta <- unique(full_umap[, .(subject_id, condition, cohort)])
donor_totals <- full_umap[, .(total_cells = .N), by = subject_id]
donor_counts <- full_umap[, .(cell_count = .N), by = .(subject_id, cell_type_analysis)]
donor_grid <- CJ(subject_id = donor_meta$subject_id, cell_type_analysis = cell_type_order, unique = TRUE)
donor_grid <- merge(donor_grid, donor_meta, by = "subject_id", all.x = TRUE, sort = FALSE)
donor_grid <- merge(donor_grid, donor_totals, by = "subject_id", all.x = TRUE, sort = FALSE)
donor_grid <- merge(donor_grid, donor_counts, by = c("subject_id", "cell_type_analysis"), all.x = TRUE, sort = FALSE)
donor_grid[is.na(cell_count), cell_count := 0L]
donor_grid[, donor_proportion := cell_count / total_cells]
donor_grid[, cell_type_label := factor(
  unname(cell_type_labels[cell_type_analysis]),
  levels = unname(cell_type_labels[cell_type_order])
)]
donor_grid[, condition := factor(as.character(condition), levels = c("HD", "SLE"))]

composition <- donor_grid[, .(
  mean_proportion = mean(donor_proportion),
  sd_proportion = sd(donor_proportion),
  se_proportion = sd(donor_proportion) / sqrt(.N),
  n_donors = .N
), by = .(condition, cell_type_analysis, cell_type_label)]
setorder(composition, condition, cell_type_label)
if (any(abs(composition[, sum(mean_proportion), by = condition]$V1 - 1) > 1e-10)) {
  stop("Donor-balanced composition does not sum to one within condition")
}

fwrite(donor_grid, file.path(record_root, "Figure5C_donor_level_cell_composition.tsv"), sep = "\t")
fwrite(composition, file.path(record_root, "Figure5C_condition_mean_cell_composition.tsv"), sep = "\t")

p5c <- ggplot(composition, aes(condition, mean_proportion)) +
  ggalluvial::geom_flow(
    aes(stratum = cell_type_label, alluvium = cell_type_label, fill = cell_type_label),
    width = 0.62, alpha = 0.28, knot.pos = 0
  ) +
  geom_col(aes(fill = cell_type_label), width = 0.62, colour = "white", linewidth = 0.22) +
  scale_fill_manual(
    values = setNames(cell_palette, unname(cell_type_labels[names(cell_palette)])),
    breaks = unname(cell_type_labels[cell_type_order]), drop = FALSE, name = "Cell type"
  ) +
  scale_x_discrete(labels = c(HD = "HD\n(n = 16 donors)", SLE = "SLE\n(n = 40 donors)")) +
  scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(expand = FALSE) +
  labs(title = "Donor-balanced cell composition", x = NULL, y = "Mean within-donor cell proportion") +
  theme_article(9.2) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "#30343A", linewidth = 0.45),
    legend.position = "right",
    legend.key.height = grid::unit(0.72, "lines"),
    legend.key.width = grid::unit(0.72, "lines"),
    legend.text = element_text(size = 7.6)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))
noop_panel(
  p5c, "C", "Figure5C_donor_balanced_cell_composition_flow",
  "Donor-balanced cell composition",
  c(full_umap_path, "reference: 14.scRNAtoolVis/14.19_ordered_barplot_by_Group.pdf"),
  7.4, 5.6
)

donor_effects_path <- "results/advanced_program_landscape/01_donor_program_atlas/donor_program_effects.tsv"
donor_effects <- read_project_tsv(donor_effects_path)
unified_ifn_effects_path <- "results/reviewer_response/04_cd14_nk_unified/all_celltype_unified_ifn_effects.tsv"
unified_ifn_effects <- read_project_tsv(unified_ifn_effects_path)
ifn_rows <- donor_effects[program == "MAGENTA_IFN_ANTIVIRAL"]
ifn_match <- match(ifn_rows$cell_type, unified_ifn_effects$cell_type)
stopifnot(!anyNA(ifn_match))
donor_effects[program == "MAGENTA_IFN_ANTIVIRAL", `:=`(
  estimate = unified_ifn_effects$estimate[ifn_match],
  std_error = unified_ifn_effects$standard_error[ifn_match],
  conf_low = unified_ifn_effects$conf_low[ifn_match],
  conf_high = unified_ifn_effects$conf_high[ifn_match],
  p_value = unified_ifn_effects$p_value[ifn_match],
  fdr = unified_ifn_effects$fdr_cell_effects[ifn_match],
  fdr_global = unified_ifn_effects$fdr_cell_effects[ifn_match],
  significant = unified_ifn_effects$fdr_cell_effects[ifn_match] < 0.05,
  model_formula = "outcome_z_within_celltype ~ condition * cell_type + age + sex + batch + cell_proportion + (1 | donor)",
  effect_scale = "within-cell-type standardized IFN score"
)]
program_id_labels <- c(
  CYAN_NEUTROPHIL_IMMATURE = "Neutrophil granule/antimicrobial",
  GREENYELLOW_B_CELL_COMPLEMENT = "B-cell/BCR",
  MAGENTA_IFN_ANTIVIRAL = "IFN/antiviral"
)
donor_effects[, program_label := factor(
  unname(program_id_labels[program]),
  levels = c("IFN/antiviral", "B-cell/BCR", "Neutrophil granule/antimicrobial")
)]
donor_effects[, cell_label := factor(
  unname(cell_type_labels[cell_type]),
  levels = rev(unname(cell_type_labels[unique(cell_type)]))
)]
effect_limit <- max(abs(donor_effects$estimate), na.rm = TRUE)
p5d <- ggplot(donor_effects, aes(program_label, cell_label, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.85) +
  geom_tile(data = donor_effects[significant == TRUE], fill = NA, colour = "#111111", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%+.2f", estimate)), size = 2.95) +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "#F7F7F5", high = "#C53030", midpoint = 0,
    limits = c(-effect_limit, effect_limit), oob = squish, name = "Adjusted\nSLE effect"
  ) +
  labs(title = "Donor-level program effects across cell types", x = NULL, y = NULL) +
  theme_article(9.3) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")
noop_panel(
  p5d, "D", "Figure5D_donor_program_effect_atlas",
  "Donor-level program effects across cell types", c(donor_effects_path, unified_ifn_effects_path), 8.4, 6.0
)

factor_scores_path <- "results/advanced_program_landscape/02_multicellular_tensor/donor_factor_scores.tsv"
factor_scores <- read_project_tsv(factor_scores_path)
factor_scores[, condition := factor(condition, levels = c("HD", "SLE"))]
p5f <- ggplot(factor_scores, aes(condition, Factor1, colour = condition, fill = condition)) +
  geom_violin(width = 0.78, trim = TRUE, alpha = 0.14, linewidth = 0) +
  geom_boxplot(width = 0.23, outlier.shape = NA, fill = "white", colour = "#34383E", linewidth = 0.55) +
  geom_point(position = position_jitter(width = 0.09, seed = 20260810), size = 2.1, alpha = 0.9) +
  scale_colour_manual(values = condition_colors, guide = "none") +
  scale_fill_manual(values = condition_colors, guide = "none") +
  labs(title = "IFN-aligned multicellular factor", x = NULL, y = "Standardized Factor 1 score") +
  theme_article(9.5)
noop_panel(
  p5f, "F", "Figure5F_IFN_aligned_factor_scores",
  "IFN-aligned multicellular factor", factor_scores_path, 5.4, 4.4
)

factor_alignment_path <- "results/advanced_program_landscape/02_multicellular_tensor/factor_program_alignment.tsv"
factor_alignment <- read_project_tsv(factor_alignment_path)
factor_alignment[, factor := factor(factor, levels = rev(paste0("Factor", 1:10)))]
factor_alignment[, program_label := factor(
  unname(program_id_labels[program]),
  levels = c("IFN/antiviral", "B-cell/BCR", "Neutrophil granule/antimicrobial")
)]
p5g <- ggplot(factor_alignment, aes(program_label, factor, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = factor_alignment[as.character(factor) == "Factor1"], fill = NA, colour = "#111111", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2.7) +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "white", high = "#C53030", midpoint = 0,
    limits = c(-1, 1), name = "Spearman rho"
  ) +
  labs(title = "Multicellular factor-program alignment", x = NULL, y = NULL) +
  theme_article(9.0) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")
noop_panel(
  p5g, "G", "Figure5G_factor_program_alignment",
  "Multicellular factor-program alignment", factor_alignment_path, 7.2, 5.4
)

factor_loadings_path <- "results/advanced_program_landscape/02_multicellular_tensor/selected_ifn_factor_loadings.tsv"
factor_loadings <- read_project_tsv(factor_loadings_path)
gene_rank <- factor_loadings[, .(
  mean_loading = mean(loading, na.rm = TRUE),
  ifn = any(is_ifn_gene)
), by = gene][order(-mean_loading)]
loading_genes <- unique(c(
  gene_rank[ifn == TRUE, head(gene, 13)],
  intersect(c("STAT1", "STAT2", "IRF7", "JAK2"), gene_rank$gene)
))
loading_data <- factor_loadings[gene %chin% loading_genes]
loading_data[, gene := factor(gene, levels = rev(gene_rank[gene %chin% loading_genes, gene]))]
loading_data[, cell_label := factor(
  unname(cell_type_labels[cell_type]),
  levels = unname(cell_type_labels[unique(cell_type)])
)]
loading_limit <- quantile(abs(loading_data$loading), 0.98, na.rm = TRUE)
p5h <- ggplot(loading_data, aes(cell_label, gene, fill = loading)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "#F7F7F5", high = "#C53030", midpoint = 0,
    limits = c(-loading_limit, loading_limit), oob = squish, name = "Factor loading"
  ) +
  labs(title = "IFN-aligned factor loadings", x = NULL, y = NULL) +
  theme_article(8.8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 32, hjust = 1), legend.position = "right")
noop_panel(
  p5h, "H", "Figure5H_IFN_factor_loadings",
  "IFN-aligned factor loadings", factor_loadings_path, 8.6, 5.4
)


small <- theme(text=element_text(size=7),plot.title=element_text(size=7.5,face='bold',margin=margin(b=3)),
  axis.title=element_text(size=7),axis.text=element_text(size=6.5),
  legend.text=element_text(size=6.5),legend.title=element_text(size=6.5),
  strip.text=element_text(size=7),plot.margin=margin(3,3,3,3),
  legend.key.size=grid::unit(2.3,'mm'),legend.spacing.x=grid::unit(1,'mm'))
shrink_text <- function(p,size=2.1) {
  for(i in seq_along(p$layers)) if(inherits(p$layers[[i]]$geom,'GeomText')) p$layers[[i]]$aes_params$size<-size
  p + small
}
p5a <- p5a + labs(title='PBMC atlas') + small + theme(legend.position='none')
atlas_labels <- full_umap[,.(label_x=median(UMAP_1),label_y=median(UMAP_2)),by=cell_type_short]
p5a <- p5a + ggrepel::geom_label_repel(data=atlas_labels,aes(label_x,label_y,label=cell_type_short),
  inherit.aes=FALSE,colour='black',fill='white',size=1.9,label.padding=grid::unit(.02,'lines'),
  label.size=0,box.padding=.05,point.padding=.03,max.overlaps=Inf,seed=20260718,
  segment.colour='#777777',segment.size=.15,min.segment.length=0)


keep <- vapply(p5b$layers,function(z) inherits(z$geom,'GeomRaster'),logical(1))
p5b$layers <- p5b$layers[keep]
p5b <- p5b + facet_wrap(~condition_label,nrow=1,labeller=as_labeller(setNames(
  sprintf('%s (%d donors)\n%s cells',group_counts$condition,group_counts$donors,scales::comma(group_counts$cells)),
  levels(full_umap$condition_label)))) +
  labs(title='HD and SLE density',subtitle=NULL,x='UMAP 1',y=NULL) + small +
  theme(strip.text=element_text(size=6,colour='white'),axis.text=element_text(colour='white',size=6),
    axis.title=element_text(colour='white',size=6.5),plot.title=element_text(colour='white',size=7.5),
    panel.spacing=grid::unit(.6,'mm'))
p5b <- p5b + viridis::scale_fill_viridis(option='magma',direction=1,limits=c(0,density_limit),
  oob=scales::squish,breaks=c(0,density_limit),labels=scales::label_number(accuracy=.001),name='Density') +
  guides(fill=guide_colourbar(barwidth=grid::unit(15,'mm'),barheight=grid::unit(2,'mm'),title.position='left')) +
  theme(legend.position='bottom',legend.text=element_text(colour='white',size=5.8),
    legend.title=element_text(colour='white',size=6),legend.background=element_rect(fill='black',colour=NA),
    legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0))
p5b <- p5b + labs(tag='B') + theme(plot.tag=element_text(colour='white',size=9,face='bold'),
  plot.tag.position=c(.01,.995))
p5c <- p5c + small + labs(title='Cell composition',y='Mean cell proportion') + theme(legend.position='none')
p5d <- shrink_text(p5d,size=1.9) + labs(title='Adjusted program effects') +
  scale_x_discrete(labels=c('IFN','B cell','Neutrophil')) +
  scale_y_discrete(labels=function(z) unname(setNames(cell_type_short,cell_type_labels)[z])) +
  theme(axis.text.x=element_text(angle=0,hjust=.5,size=6),axis.text.y=element_text(size=6),legend.position='right') +
  guides(fill=guide_colourbar(title='Effect',barheight=grid::unit(15,'mm'),barwidth=grid::unit(2,'mm'),title.position='top'))
p5f <- p5f + small + labs(title='IFN-aligned factor',y='Standardized Factor 1')
p5g <- shrink_text(p5g,size=1.9) + labs(title='Factor-program alignment') +
  scale_x_discrete(labels=c('IFN','B cell','Neutrophil')) +
  theme(axis.text.x=element_text(angle=0,hjust=.5,size=6),axis.text.y=element_text(size=6),legend.position='right') +
  guides(fill=guide_colourbar(title='rho',barheight=grid::unit(15,'mm'),barwidth=grid::unit(2,'mm'),title.position='top'))

loading_rank <- factor_loadings[,.(mean_loading=mean(loading)),by=gene][order(-mean_loading)]
balanced_genes <- unique(c(head(loading_rank$gene,10L),tail(loading_rank$gene,10L)))
balanced <- copy(factor_loadings[gene %in% balanced_genes])
gene_order <- balanced[,.(mean_loading=mean(loading)),by=gene][order(-mean_loading)]$gene
balanced[,gene:=factor(gene,levels=rev(gene_order))]
balanced[,cell_label:=factor(unname(cell_type_labels[cell_type]),levels=unname(cell_type_labels[unique(cell_type)]))]
lim <- max(abs(balanced$loading))
p5h <- ggplot(balanced,aes(cell_label,gene,fill=loading)) +
  geom_tile(colour='white',linewidth=.15) +
  scale_fill_gradient2(low='#2B6CB0',mid='white',high='#C53030',limits=c(-lim,lim),name='Factor loading') +
  labs(title='Positive and negative Factor 1 loadings',x=NULL,y=NULL) +
  theme_article(7) + small +
  scale_x_discrete(labels=function(z) unname(setNames(cell_type_short,cell_type_labels)[z])) +
  theme(panel.grid=element_blank(),axis.text.x=element_text(angle=30,hjust=1,size=6.5),
    axis.text.y=element_text(size=5.8),legend.position='right') +
  guides(fill=guide_colourbar(barheight=grid::unit(20,'mm'),barwidth=grid::unit(2,'mm')))


layout <- 'AB\nCD\nEF\nGG'


fig <- wrap_plots(A=wrap_elements(full=p5a),B=wrap_elements(full=p5b),
  C=wrap_elements(full=p5c),D=wrap_elements(full=p5d),
  E=wrap_elements(full=p5f),F=wrap_elements(full=p5g),G=wrap_elements(full=p5h),
  design=layout,heights=c(1.65,1.30,1.30,2.25)) +
  plot_annotation(tag_levels=list(c('A','','C','D','E','F','G'))) &
  theme(plot.tag=element_text(family=font_family,size=9,face='bold'),plot.tag.position=c(0,.995))
ggsave(file.path(stage,'work/Figure5.pdf'),fig,width=5.12,height=6.5,units='in',device=cairo_pdf,family=font_family,bg='white')
ggsave(file.path(stage,'work/Figure5.png'),fig,width=5.12,height=6.5,units='in',device=ragg::agg_png,dpi=450,bg='white')
data.table::fwrite(balanced,file.path(stage,'figures/Figure5_loading_plotted_values.tsv'),sep='\t')
message('Figure5 rendered with unchanged model outputs and frozen UMAP coordinates.')
