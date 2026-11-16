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

clinical_path <- "results/exploratory_program_landscape/01_pairwise_clinical/program_high_vs_immune_low_clinical_contrasts.tsv"
clinical <- read_tsv(clinical_path)[trait %chin% c("SLEDAI", "C3", "C4", "anti_dsDNA")]
clinical[, state_label := factor(unname(state_labels[high_group]),
  levels = c("Neutrophil-high", "IFN-high", "B-cell-high"))]
clinical[, trait_label := factor(trait, levels = c("SLEDAI", "anti_dsDNA", "C3", "C4"))]
p4a <- ggplot(clinical, aes(trait_label, state_label, fill = activity_aligned_effect)) +
  geom_tile(colour = "white", linewidth = 0.9) +
  geom_tile(data = clinical[fdr < 0.05], fill = NA, colour = "#171717", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", activity_aligned_effect)), size = 3.4) +
  scale_fill_gradient(low = "#F4F5F5", high = "#B33A3A", limits = c(0, max(clinical$activity_aligned_effect)),
                      name = "Activity-aligned\neffect") +
  labs(title = "Clinical effects relative to the immune-low region", x = NULL, y = NULL) +
  theme_article(10) +
  theme(panel.grid = element_blank(), legend.position = "right")
save_plot(p4a, "Figure4", "A", "Figure4A_program_clinical_effect_map",
          "Clinical effects relative to the immune-low region", clinical_path, 8.4, 4.8)

ranking_path <- "results/exploratory_program_landscape/01_pairwise_clinical/clinical_separation_ranking.tsv"
ranking <- read_tsv(ranking_path)
ranking[, state_label := factor(unname(state_labels[high_group]), levels = rev(c("Neutrophil-high", "IFN-high", "B-cell-high")))]
p4b <- ggplot(ranking, aes(mean_activity_aligned_effect, state_label, fill = state_label)) +
  geom_col(width = 0.60) +
  geom_text(aes(label = sprintf("%.2f", mean_activity_aligned_effect)), hjust = -0.16, size = 3.6) +
  scale_fill_manual(values = state_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(title = "Clinical separation across disease traits", x = "Mean activity-aligned standardized effect", y = NULL) +
  theme_article(10)
save_plot(p4b, "Figure4", "B", "Figure4B_clinical_separation_ranking",
          "Clinical separation across disease traits", ranking_path, 7.0, 4.8)

fig4 <- p4a | p4b
combine_and_save(fig4, "Figure4", "Clinical topology of the landscape",
                 c(clinical_path, ranking_path), 14.8, 5.8)



stage <- 'F:/program/Hei/SLE_R1/checks/r1_revision_20261116/work'
dir.create(stage,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(stage,'Figure4.pdf'),fig4+plot_annotation(tag_levels='A'),width=5.12,height=4.5,units='in',device=cairo_pdf,family=font_family,bg='white')
