if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
options(stringsAsFactors = FALSE, scipen = 99)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrastr)
})

project_root <- "F:/program/Hei/SLE_v1"
args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args)) args[[1]] else "F:/program/Hei/SLE/redraw_20260810"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
for (i in 1:9) dir.create(file.path(output_root, paste0("Figure", i)), recursive = TRUE, showWarnings = FALSE)

font_family <- "Times New Roman"
program_ids <- c("magenta", "greenyellow", "cyan")
program_labels <- c(
  magenta = "IFN/antiviral",
  greenyellow = "B-cell/BCR",
  cyan = "Immature-neutrophil"
)
program_colors <- c(
  "IFN/antiviral" = "#A23B3B",
  "B-cell/BCR" = "#2F7D67",
  "Immature-neutrophil" = "#D18A27"
)
state_labels <- c(
  HEALTHY_CONTROL = "Healthy control",
  S01_NEUTROPHIL_HIGH = "Neutrophil-high",
  S02_IMMUNE_LOW = "Immune-low",
  S03_B_CELL_HIGH = "B-cell-high",
  S04_IFN_HIGH = "IFN-high"
)
state_colors <- c(
  "Healthy control" = "#3A9D78",
  "Neutrophil-high" = "#D18A27",
  "Immune-low" = "#707A86",
  "B-cell-high" = "#2F7D67",
  "IFN-high" = "#A23B3B"
)
condition_colors <- c(HD = "#2C6E9B", SLE = "#C23B3B")

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
theme_set(theme_article())
update_geom_defaults("text", list(family = font_family, colour = "#202124"))
update_geom_defaults("label", list(family = font_family, colour = "#202124"))

manifest <- list()
add_manifest <- function(figure, panel, title, sources, width, height, pdf, png, kind = "panel") {
  manifest[[length(manifest) + 1L]] <<- data.table(
    figure = figure, panel = panel, kind = kind, title = title,
    source_files = paste(sources, collapse = " | "), width_in = width, height_in = height,
    pdf = normalizePath(pdf, winslash = "/", mustWork = FALSE),
    png = normalizePath(png, winslash = "/", mustWork = FALSE)
  )
}

save_plot <- function(plot, figure, panel, stem, title, sources, width, height, kind = "panel", png_dpi = 320) {
  out_dir <- file.path(output_root, figure)
  pdf <- file.path(out_dir, paste0(stem, ".pdf"))
  png <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(
    pdf, plot, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, family = font_family, bg = "white"
  )
  ggsave(
    png, plot, width = width, height = height, units = "in", dpi = png_dpi,
    device = ragg::agg_png, background = "white"
  )
  add_manifest(figure, panel, title, sources, width, height, pdf, png, kind)
  invisible(c(pdf = pdf, png = png))
}

combine_and_save <- function(plot, figure, title, sources, width, height) {
  tagged <- plot + plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(family = font_family, face = "bold", size = 15),
      plot.tag.position = c(0.002, 0.998)
    )
  save_plot(
    tagged, figure, "combined", paste0(figure, "_combined"), title, sources,
    width, height, kind = "combined"
  )
  invisible(tagged)
}

read_tsv <- function(...) fread(file.path(project_root, ...), sep = "\t", na.strings = c("", "NA"))
fmt_p <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x))

cell_type_labels <- c(
  B = "B", B_Naive = "Naive B", Plasmablast = "Plasmablast", Monocyte = "Monocyte",
  CD14_Monocyte = "CD14 monocyte", CD16_Monocyte = "CD16 monocyte", pDC = "pDC",
  CD4_T = "CD4 T", CD4_T_Naive = "Naive CD4 T", CD8_T_Naive = "Naive CD8 T",
  CD8_T_GZMK = "GZMK CD8 T", CD8_T_Cytotoxic = "Cytotoxic CD8 T", NK = "NK",
  Platelet = "Platelet", Erythroid = "Erythroid"
)

hubs_path <- "results/3.1/06_wgcna/hub_genes.tsv"
bcell_summary_path <- "results/exploratory_program_landscape/04_followup_mechanisms/01_bcell_localization/bcell_program_celltype_summary.tsv"
bcell_summary <- read_tsv(bcell_summary_path)
bcell_localization <- bcell_summary[, .(
  mean_score = weighted.mean(mean_score, pmax(n_profiles, 1), na.rm = TRUE)
), by = .(cohort, cell_type_analysis)]
bcell_localization[, cell_label := unname(cell_type_labels[cell_type_analysis])]
bcell_localization[is.na(cell_label), cell_label := cell_type_analysis]
cell_rank <- bcell_localization[, .(score = mean(mean_score)), by = cell_label][order(score), cell_label]
bcell_localization[, cell_label := factor(cell_label, levels = cell_rank)]
p8a <- ggplot(bcell_localization, aes(mean_score, cell_label, colour = cohort)) +
  geom_segment(aes(x = 0, xend = mean_score, yend = cell_label), colour = "#D8DADD", linewidth = 0.48,
               position = position_dodge(width = 0.48)) +
  geom_point(position = position_dodge(width = 0.48), size = 2.7) +
  scale_colour_manual(values = c(Child = "#4C78A8", Adult = "#E45756"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(title = "B-cell program localization across cell types", x = "Mean B-cell/complement program score", y = NULL) +
  theme_article(9.3)
save_plot(p8a, "Figure8", "A", "Figure8A_Bcell_program_localization",
          "B-cell program localization across cell types", bcell_summary_path, 8.2, 6.2)

bcell_decomp_path <- "results/exploratory_program_landscape/04_followup_mechanisms/01_bcell_localization/bcell_gene_set_decomposition.tsv"
bcell_decomp <- read_tsv(bcell_decomp_path)
identity_genes_text <- bcell_decomp[gene_set == "B_cell_identity", genes][1]
identity_genes <- unique(strsplit(identity_genes_text, ";", fixed = TRUE)[[1]])
greenyellow_hubs <- read_tsv(hubs_path)[module == "greenyellow"][order(hub_rank)][1:18]
greenyellow_hubs[, category := fifelse(gene %chin% identity_genes, "B-cell identity gene", "Other program gene")]
greenyellow_hubs[, gene := factor(gene, levels = gene[order(kME)])]
p8b <- ggplot(greenyellow_hubs, aes(kME, gene, colour = category)) +
  geom_segment(aes(x = 0, xend = kME, yend = gene), colour = "#D8DADD", linewidth = 0.5) +
  geom_point(size = 2.7) +
  scale_colour_manual(values = c("B-cell identity gene" = "#2F7D67", "Other program gene" = "#A7ADB4"), name = NULL) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "B-cell program hub genes", x = "Module membership (kME)", y = NULL) +
  theme_article(9.0)
save_plot(p8b, "Figure8", "B", "Figure8B_Bcell_program_hub_genes",
          "B-cell program hub genes", c(hubs_path, bcell_decomp_path), 6.8, 6.2)

decomp_rows <- bcell_decomp[cell_label == "B_Naive"]
decomp_wide <- dcast(decomp_rows, comparison ~ gene_set, value.var = "n_genes_present")
decomp_plot <- decomp_wide[, .(
  comparison,
  `B-cell identity` = B_cell_identity,
  `BAFF/NF-kB associated` = BAFF_NFkB_associated,
  `Other program genes` = full_B_cell_complement_program - B_cell_identity - BAFF_NFkB_associated
)]
decomp_plot <- melt(decomp_plot, id.vars = "comparison", variable.name = "category", value.name = "genes")
decomp_plot[, cohort := factor(ifelse(grepl("^Child", comparison), "Child", "Adult"), levels = c("Child", "Adult"))]
decomp_plot[, category := factor(category, levels = c("Other program genes", "B-cell identity", "BAFF/NF-kB associated"))]
p8c <- ggplot(decomp_plot, aes(cohort, genes, fill = category)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = genes), position = position_stack(vjust = 0.5), size = 3.3) +
  scale_fill_manual(values = c("Other program genes" = "#B8BEC5", "B-cell identity" = "#2F7D67", "BAFF/NF-kB associated" = "#D18A27"), name = NULL) +
  labs(title = "Detected program-gene composition", x = NULL, y = "Detected genes") +
  theme_article(9.3)
save_plot(p8c, "Figure8", "C", "Figure8C_Bcell_program_gene_composition",
          "Detected program-gene composition", bcell_decomp_path, 6.8, 4.6)

fig8 <- (p8a | ((p8b / p8c) + plot_layout(heights = c(1.2, 0.8)))) +
  plot_layout(widths = c(1.08, 0.92))
combine_and_save(fig8, "Figure8", "B-cell identity state",
                 c(bcell_summary_path, hubs_path, bcell_decomp_path), 15.2, 10.0)



stage <- 'F:/program/Hei/SLE_R1/checks/r1_revision_20261116/work'
dir.create(stage,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(stage,'Figure7.pdf'),fig8+plot_annotation(tag_levels='A'),width=5.12,height=4.5,units='in',device=cairo_pdf,family=font_family,bg='white')
