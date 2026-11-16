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
output_root <- if (length(args)) args[[1]] else "F:/program/Hei/SLE/figure2_maintext_stage/output"
figure_dir <- file.path(output_root, "Figure2")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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
add_manifest <- function(panel, title, sources, width, height, pdf, png, kind = "panel") {
  manifest[[length(manifest) + 1L]] <<- data.table(
    figure = "Figure2", panel = panel, kind = kind, title = title,
    source_files = paste(sources, collapse = " | "), width_in = width, height_in = height,
    pdf = normalizePath(pdf, winslash = "/", mustWork = FALSE),
    png = normalizePath(png, winslash = "/", mustWork = FALSE)
  )
}

save_plot <- function(plot, panel, stem, title, sources, width, height, kind = "panel", png_dpi = 320) {
  pdf <- file.path(figure_dir, paste0(stem, ".pdf"))
  png <- file.path(figure_dir, paste0(stem, ".png"))
  ggsave(pdf, plot, width = width, height = height, units = "in",
         device = grDevices::cairo_pdf, family = font_family, bg = "white")
  ggsave(png, plot, width = width, height = height, units = "in", dpi = png_dpi,
         device = ragg::agg_png, background = "white")
  add_manifest(panel, title, sources, width, height, pdf, png, kind)
}

read_tsv <- function(path) fread(file.path(project_root, path), sep = "\t", na.strings = c("", "NA"))
fmt_p <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x))

deg_path <- "results/3.2/02_deg_gsea/full_transcriptome_limma.tsv"
gsea_path <- "results/3.2/02_deg_gsea/fgsea_all_collections.tsv"
pres_path <- "results/3.2/03_preservation/module_preservation_statistics.tsv"
concord_path <- "results/3.2/03_preservation/module_membership_concordance.tsv"


deg <- read_tsv(deg_path)
deg[, fdr_plot := -log10(pmax(fdr, 1e-300))]
deg[, direction_plot := factor(direction, levels = c("SLE_up", "SLE_down", "not_primary"))]
deg_counts <- deg[, .N, by = direction]
up_n <- deg_counts[direction == "SLE_up", N]
down_n <- deg_counts[direction == "SLE_down", N]
label_data <- rbindlist(list(
  deg[direction == "SLE_up"][order(-logFC)][1:10],
  deg[direction == "SLE_down"][order(logFC)][1:6]
), use.names = TRUE)

p2a <- ggplot(deg, aes(logFC, fdr_plot, colour = direction_plot)) +
  ggrastr::rasterise(geom_point(alpha = 0.38, size = 0.80), dpi = 320) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted", colour = "#666666", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "#8A8A8A", linewidth = 0.45) +
  ggrepel::geom_text_repel(
    data = label_data, aes(label = gene), family = font_family,
    size = 3.1, colour = "#222222", box.padding = 0.35, point.padding = 0.2,
    min.segment.length = 0, max.overlaps = Inf, seed = 20260810, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(SLE_up = "#C73E1D", SLE_down = "#2166AC", not_primary = "#B8B8B8"),
    breaks = c("SLE_up", "SLE_down", "not_primary"),
    labels = c(paste0("SLE-up (n=", up_n, ")"), paste0("SLE-down (n=", down_n, ")"), "Other"),
    name = NULL
  ) +
  labs(
    title = "Adult SLE differential-expression landscape",
    x = expression(log[2]~"fold change (SLE vs healthy)"),
    y = expression(-log[10]~"FDR")
  ) +
  theme_article(10) +
  theme(panel.grid.major.x = element_blank())
save_plot(p2a, "A", "Figure2A_adult_full_transcriptome_volcano",
          "Adult SLE differential-expression landscape", deg_path, 7.5, 5.5)


gsea_all <- read_tsv(gsea_path)
gsea_landscape <- gsea_all[is.finite(padj) & (NES > 0 | collection == "FROZEN")][
  order(padj, -NES), head(.SD, if (collection[[1]] == "FROZEN") 3L else 14L), by = collection
]

clean_gsea_label <- function(label, collection) {
  if (collection == "FROZEN") {
    return(switch(label,
      MAGENTA_IFN_ANTIVIRAL = "IFN/antiviral program",
      GREENYELLOW_B_CELL_COMPLEMENT = "B-cell/BCR program",
      CYAN_NEUTROPHIL_IMMATURE = "Immature-neutrophil program",
      label
    ))
  }
  ontology <- if (grepl("^GOBP_", label)) {
    "GO BP"
  } else if (grepl("^GOCC_", label)) {
    "GO CC"
  } else if (grepl("^GOMF_", label)) {
    "GO MF"
  } else {
    "KEGG"
  }
  core <- sub("^(GOBP_|GOCC_|GOMF_|KEGG_)", "", label)
  core <- tools::toTitleCase(tolower(gsub("_", " ", core)))
  paste0(core, " [", ontology, "]")
}

gsea_landscape[, collection_label := factor(collection,
  levels = c("FROZEN", "GO", "KEGG"), labels = c("Frozen programs", "Gene Ontology", "KEGG"))]
gsea_landscape[, display_label := mapply(clean_gsea_label, label, collection, USE.NAMES = FALSE)]
gsea_landscape[, plot_key := paste(collection, display_label, sep = "___")]
gsea_landscape <- gsea_landscape[order(collection_label, NES)]
gsea_landscape[, plot_key := factor(plot_key, levels = unique(plot_key))]
gsea_landscape[, fdr_strength := pmin(-log10(pmax(padj, 1e-100)), 100)]

p2b <- ggplot(gsea_landscape, aes(NES, plot_key, size = fdr_strength, fill = NES)) +
  geom_vline(xintercept = 0, colour = "#AEB3B8", linewidth = 0.5) +
  geom_point(shape = 21, alpha = 0.95, colour = "#4B5157", stroke = 0.25) +
  facet_grid(collection_label ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, name = "NES") +
  scale_size_continuous(range = c(2.2, 6.3), name = expression(-log[10]~"FDR")) +
  labs(title = "Adult preranked GSEA landscape", x = "Normalized enrichment score", y = NULL) +
  theme_article(9.2) +
  theme(
    panel.grid.major.x = element_line(colour = "#ECECEC", linewidth = 0.3),
    legend.box = "horizontal",
    legend.position = "bottom"
  )
save_plot(p2b, "B", "Figure2B_adult_GSEA_landscape",
          "Adult preranked GSEA landscape", gsea_path, 13.2, 10.4)


program_key <- c(
  MAGENTA_IFN_ANTIVIRAL = "magenta",
  GREENYELLOW_B_CELL_COMPLEMENT = "greenyellow",
  CYAN_NEUTROPHIL_IMMATURE = "cyan"
)
gsea_program <- gsea_all[collection == "FROZEN" & NES > 0 & padj < 0.05]
gsea_program[, program := factor(unname(program_labels[program_key[label]]),
                                 levels = rev(unname(program_labels)))]
p2c <- ggplot(gsea_program, aes(NES, program, fill = program)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = paste0("NES=", sprintf("%.2f", NES), "  FDR=", fmt_p(padj))),
            hjust = -0.05, size = 3.2) +
  scale_fill_manual(values = program_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.42))) +
  labs(title = "Adult enrichment of pediatric programs", x = "Normalized enrichment score", y = NULL) +
  theme_article(10)
save_plot(p2c, "C", "Figure2C_adult_program_enrichment",
          "Adult enrichment of pediatric programs", gsea_path, 7.5, 5.5)


pres <- read_tsv(pres_path)
pres[, selected := module %chin% program_ids]
pres[, module := factor(module, levels = pres[order(Zsummary.pres), module])]
p2d <- ggplot(pres, aes(Zsummary.pres, module, colour = selected)) +
  geom_vline(xintercept = c(2, 10), linetype = c("dotted", "dashed"), colour = "#9AA0A6", linewidth = 0.55) +
  geom_segment(aes(x = 0, xend = Zsummary.pres, yend = module), colour = "#D8DADD", linewidth = 0.5) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = c(`FALSE` = "#8A939D", `TRUE` = "#A23B3B"), guide = "none") +
  labs(title = "Module preservation in adult SLE", x = "Preservation Zsummary", y = NULL) +
  theme_article(9.5)
save_plot(p2d, "D", "Figure2D_adult_module_preservation",
          "Module preservation in adult SLE", pres_path, 7.2, 5.5)


concord <- read_tsv(concord_path)
concord[, selected := module %chin% program_ids]
concord[, module := factor(module, levels = concord[order(spearman_kME), module])]
p2e <- ggplot(concord, aes(spearman_kME, module, colour = selected)) +
  geom_vline(xintercept = 0, colour = "#B5B8BC", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = spearman_kME, yend = module), colour = "#D8DADD", linewidth = 0.5) +
  geom_point(aes(size = n_genes), alpha = 0.95) +
  scale_colour_manual(values = c(`FALSE` = "#8A939D", `TRUE` = "#2F7D67"), guide = "none") +
  scale_size_continuous(range = c(2.2, 5.5), name = "Genes") +
  labs(title = "Cross-age module-membership concordance", x = "Spearman correlation of kME", y = NULL) +
  theme_article(9.5)
save_plot(p2e, "E", "Figure2E_cross_age_kME_concordance",
          "Cross-age module-membership concordance", concord_path, 7.2, 5.5)

fig2 <- p2a / p2b / p2c / (p2d | p2e) +
  plot_layout(heights = c(1.05, 1.80, 0.72, 1.05))
tagged <- fig2 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = font_family, face = "bold", size = 15),
        plot.tag.position = c(0.002, 0.998))
save_plot(tagged, "combined", "Figure2_combined", "Adult differential landscape and cross-age validation",
          c(deg_path, gsea_path, pres_path, concord_path), 15.8, 22.0, kind = "combined", png_dpi = 220)

fwrite(rbindlist(manifest, fill = TRUE), file.path(output_root, "figure2_manifest.tsv"), sep = "\t", quote = TRUE)

cat(sprintf("DEG_COUNTS\tall=%d\tSLE_up=%d\tSLE_down=%d\n", nrow(deg), up_n, down_n))
cat(sprintf("GSEA_ROWS\tselected=%d\tFROZEN=%d\tGO=%d\tKEGG=%d\n",
            nrow(gsea_landscape),
            gsea_landscape[collection == "FROZEN", .N],
            gsea_landscape[collection == "GO", .N],
            gsea_landscape[collection == "KEGG", .N]))

stage <- 'F:/program/Hei/SLE_R1/checks/r1_revision_20261116/work'
dir.create(stage,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(stage,'Figure2.pdf'),tagged,width=5.12,height=7.4,units='in',device=cairo_pdf,family=font_family,bg='white')
