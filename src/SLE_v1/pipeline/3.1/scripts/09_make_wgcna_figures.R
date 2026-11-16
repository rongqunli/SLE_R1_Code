if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "3.1")
figure_dir <- file.path(results_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The data.table and ggplot2 packages are required for WGCNA figures.")
}

data.table::setDTthreads(1L)
read_tsv <- function(path) data.table::fread(path, sep = "\t", na.strings = "NA")
save_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
}
theme_wgcna <- function(size = 10) {
  ggplot2::theme_classic(base_size = size) + ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(color = "#444444"),
    axis.text = ggplot2::element_text(color = "#252525"),
    panel.grid = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )
}

program_labels <- c(magenta = "IFN/antiviral", greenyellow = "B-cell/complement", cyan = "Immature-neutrophil")
program_colors <- c("IFN/antiviral" = "#A23B3B", "B-cell/complement" = "#2F7D67", "Immature-neutrophil" = "#D18A27")
all_traits <- read_tsv(file.path(results_root, "06_wgcna", "module_trait_associations.tsv"))
traits <- data.table::copy(all_traits)
traits <- traits[model_family %in% c("activity_composition_adjusted", "renal_composition_adjusted", "c3_composition_adjusted", "c4_composition_adjusted")]
traits[, aligned_t := ifelse(grepl("^c[34]_", model_family), -t_statistic, t_statistic)]
traits[, axis := c(activity_composition_adjusted = "SLEDAI adjusted", renal_composition_adjusted = "Renal activity adjusted", c3_composition_adjusted = "Low C3 adjusted", c4_composition_adjusted = "Low C4 adjusted")[model_family]]
audit <- read_tsv(file.path(results_root, "07_module_selection_audit", "Table_S2_module_selection_audit.tsv"))
stability <- read_tsv(file.path(results_root, "06_wgcna", "module_stability.tsv"))

selection <- traits[, .(selection_score = max(aligned_t), best_axis = axis[which.max(aligned_t)], minimum_adjusted_fdr = min(fdr), significant_adjusted_axes = sum(fdr < 0.05)), by = module]
selection <- merge(selection, stability, by = "module", all.x = TRUE)
selection[, eligible := significant_adjusted_axes >= 1 & median_kme_rank_correlation >= 0.97 & proportion_correlation_ge_0_7 == 1]
selection[, positive_candidate := eligible & selection_score > 0]
data.table::setorder(selection, -positive_candidate, -selection_score, module)
selection[, selection_rank := NA_integer_]
selection[positive_candidate == TRUE, selection_rank := seq_len(.N)]
selection[, selected := positive_candidate & selection_rank <= 3]
selection[, biological_label := unname(program_labels[module])]
traits <- merge(traits, selection[, .(module, selected, selection_rank, biological_label)], by = "module", all.x = TRUE)
module_order <- selection[order(selection_score), module]
traits[, module_display := ifelse(selected, paste0(selection_rank, ". ", module, " | ", biological_label), module)]
display_levels <- ifelse(selection$selected[match(module_order, selection$module)], paste0(selection$selection_rank[match(module_order, selection$module)], ". ", module_order, " | ", selection$biological_label[match(module_order, selection$module)]), module_order)
traits[, module_display := factor(module_display, levels = display_levels)]
traits[, cell_label := paste0("t=", sprintf("%+.1f", aligned_t), "\nq=", ifelse(fdr < 0.001, format(fdr, scientific = TRUE, digits = 1), sprintf("%.3f", fdr)))]

p_heatmap <- ggplot2::ggplot(traits, ggplot2::aes(x = axis, y = module_display, fill = pmax(-8, pmin(8, aligned_t)))) +
  ggplot2::geom_tile(color = "white", linewidth = 0.55) +
  ggplot2::geom_tile(data = traits[selected == TRUE], fill = NA, color = "#202020", linewidth = 0.9) +
  ggplot2::geom_text(ggplot2::aes(label = cell_label), size = 3.1, lineheight = 0.92) +
  ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C84343", midpoint = 0, limits = c(-8, 8), name = "Activity-aligned\nt statistic") +
  ggplot2::labs(title = "Adjusted module-trait landscape", subtitle = "Outlined rows are the three selected transport programs", x = NULL, y = NULL) +
  theme_wgcna(10) + ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", angle = 25, hjust = 1), axis.text.y = ggplot2::element_text(size = 8.5))
save_plot(p_heatmap, "Figure1C_module_selection_heatmap", 12, 6.2)

full_heatmap <- data.table::copy(all_traits)
full_heatmap[, aligned_t := ifelse(grepl("^c[34]_", model_family), -t_statistic, t_statistic)]
family_order <- c("activity_base", "activity_composition_adjusted", "renal_base", "renal_composition_adjusted", "c3_base", "c3_composition_adjusted", "c4_base", "c4_composition_adjusted")
family_labels <- c(
  activity_base = "Activity base", activity_composition_adjusted = "Activity + composition",
  renal_base = "Renal base", renal_composition_adjusted = "Renal + composition",
  c3_base = "C3 base", c3_composition_adjusted = "C3 + composition",
  c4_base = "C4 base", c4_composition_adjusted = "C4 + composition"
)
full_heatmap[, family_display := factor(unname(family_labels[model_family]), levels = unname(family_labels[family_order]))]
full_heatmap[, significance := ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", "")))]
full_heatmap[, module := factor(module, levels = rev(selection$module))]
p_full_heatmap <- ggplot2::ggplot(full_heatmap, ggplot2::aes(x = family_display, y = module, fill = pmax(-8, pmin(8, aligned_t)))) +
  ggplot2::geom_tile(color = "white", linewidth = 0.45) + ggplot2::geom_text(ggplot2::aes(label = significance), size = 3.4) +
  ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C84343", midpoint = 0, limits = c(-8, 8), name = "Activity-aligned\nt statistic") +
  ggplot2::labs(title = "WGCNA modules track distinct clinical axes", subtitle = "Positive values indicate higher module expression with greater disease activity; stars indicate BH FDR thresholds", x = NULL, y = "Module") +
  theme_wgcna(10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), panel.grid = ggplot2::element_blank())
save_plot(p_full_heatmap, "FigureS_module_trait_heatmap", 12, 8)

stability_plot <- stability[module %in% names(program_labels)]
stability_plot[, program := factor(unname(program_labels[module]), levels = rev(unname(program_labels)))]
p_stability <- ggplot2::ggplot(stability_plot, ggplot2::aes(x = median_kme_rank_correlation, y = program, color = program)) +
  ggplot2::geom_segment(ggplot2::aes(x = q025_kme_rank_correlation, xend = q975_kme_rank_correlation, yend = program), linewidth = 1.2) +
  ggplot2::geom_point(size = 3.2) + ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n_genes), x = q975_kme_rank_correlation + 0.0008), hjust = 0, color = "#252525", size = 3.1) +
  ggplot2::scale_color_manual(values = program_colors, guide = "none") + ggplot2::scale_x_continuous(limits = c(0.968, 0.997), breaks = seq(0.97, 0.995, 0.005)) +
  ggplot2::labs(title = "Module membership stability across patient bootstraps", subtitle = "Median kME rank correlation; bars show 2.5th-97.5th percentiles across 25 bootstraps", x = "Patient-bootstrap kME rank correlation", y = NULL) + theme_wgcna(10)
save_plot(p_stability, "Figure1D_module_stability", 7.5, 4.2)

hubs <- read_tsv(file.path(results_root, "06_wgcna", "hub_genes.tsv"))
hubs <- hubs[module %in% names(program_labels)][order(module, hub_rank)][, head(.SD, 8), by = module]
hubs[, program := factor(unname(program_labels[module]), levels = unname(program_labels))]
hubs[, gene_panel := paste(gene, program, sep = "___")]
hubs[, gene_panel := factor(gene_panel, levels = gene_panel[order(kME)]), by = program]
p_hubs <- ggplot2::ggplot(hubs, ggplot2::aes(x = kME, y = gene_panel, color = program)) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = kME, yend = gene_panel), color = "#D5D8DC", linewidth = 0.5) + ggplot2::geom_point(size = 2.7) +
  ggplot2::facet_wrap(~program, nrow = 1, scales = "free_y") + ggplot2::scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  ggplot2::scale_color_manual(values = program_colors, guide = "none") + ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  ggplot2::labs(title = "Representative hub genes", subtitle = "Top eight genes by within-module kME rank", x = "Module membership (kME)", y = NULL) + theme_wgcna(9)
save_plot(p_hubs, "Figure1E_representative_hubs", 10, 5.2)

data.table::fwrite(data.table::data.table(
  figure = c("Figure1C", "Figure1D", "Figure1E", "FigureS"),
  file_png = c("Figure1C_module_selection_heatmap.png", "Figure1D_module_stability.png", "Figure1E_representative_hubs.png", "FigureS_module_trait_heatmap.png"),
  file_pdf = c("Figure1C_module_selection_heatmap.pdf", "Figure1D_module_stability.pdf", "Figure1E_representative_hubs.pdf", "FigureS_module_trait_heatmap.pdf"),
  source = c("module_trait_associations.tsv + module_selection_audit", "module_stability.tsv", "hub_genes.tsv", "module_trait_associations.tsv")
), file.path(figure_dir, "figure_manifest.tsv"), sep = "\t", quote = FALSE)
message("3.1 WGCNA figures written to ", figure_dir)
