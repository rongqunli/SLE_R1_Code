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
})

project_root <- "F:/program/Hei/SLE_v1"
args <- c(file.path(stage,'figures/Figure8_panels'))
output_dir <- if (length(args)) args[[1]] else file.path(project_root, "Figures", "Figure9")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Times New Roman"
theme_article <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, colour = "#202124"),
      plot.title = element_text(face = "bold", size = rel(1.16), margin = margin(b = 7)),
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

read_tsv <- function(relative_path) {
  fread(file.path(project_root, relative_path), sep = "\t", na.strings = c("", "NA"))
}

save_panel <- function(plot, stem, width, height) {
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, family = font_family, bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")), plot,
    width = width, height = height, units = "in", dpi = 320,
    device = ragg::agg_png, background = "white"
  )
}

flow_path <- "results/reviewer_response/08_trial_reporting_ancova_permutation/trial_flow_parent_and_molecular.tsv"
ancova_components_path <- "results/reviewer_response/08_trial_reporting_ancova_permutation/ancova_coordinate_effects_trial_specific.tsv"
moderation_trial_path <- "results/reviewer_response/08_trial_reporting_ancova_permutation/ancova_bcell_baseline_interactions_trial_specific.tsv"
moderation_meta_path <- "results/reviewer_response/08_trial_reporting_ancova_permutation/ancova_bcell_baseline_interactions_fixed_effect_meta.tsv"
abundance_samples_path <- "results/reviewer_response/07_bcell_abundance/sample_level_b_abundance_and_program_scores.tsv"
abundance_correlations_path <- "results/reviewer_response/07_bcell_abundance/b_program_abundance_correlations.tsv"
abundance_models_path <- "results/reviewer_response/07_bcell_abundance/abundance_adjusted_treatment_models.tsv"

ancova_components <- read_tsv(ancova_components_path)
coord_names <- c(
  "coord_CYAN_NEUTROPHIL_IMMATURE",
  "coord_GREENYELLOW_B_CELL_COMPLEMENT",
  "coord_MAGENTA_IFN_ANTIVIRAL"
)
component_labels <- c(
  coord_CYAN_NEUTROPHIL_IMMATURE = "Immature-neutrophil",
  coord_GREENYELLOW_B_CELL_COMPLEMENT = "B-cell-associated",
  coord_MAGENTA_IFN_ANTIVIRAL = "IFN/antiviral"
)
dose_labels <- c(Tabalumab_Q2W = "Q2W", Tabalumab_Q4W = "Q4W")
trial_colors <- c("ILLUMINATE-1" = "#445E93", "ILLUMINATE-2" = "#C4515C")

flow <- read_tsv(flow_path)
flow_base <- unique(flow[, .(
  trial, treatment, parent_trial_randomized,
  transcriptomic_cohort_subjects, baseline_available
)])
flow_long <- rbind(
  flow_base[, .(trial, treatment, parent_trial_randomized, stage = "parent_trial_randomized", participants = parent_trial_randomized)],
  flow_base[, .(trial, treatment, parent_trial_randomized, stage = "transcriptomic_cohort_subjects", participants = transcriptomic_cohort_subjects)],
  flow_base[, .(trial, treatment, parent_trial_randomized, stage = "baseline_available", participants = baseline_available)],
  flow[followup == "week16", .(
    trial, treatment, parent_trial_randomized,
    stage = "paired_week16", participants = paired_baseline_followup
  )],
  flow[followup == "week52", .(
    trial, treatment, parent_trial_randomized,
    stage = "paired_week52", participants = paired_baseline_followup
  )]
)
flow_long[, stage := factor(
  stage,
  levels = c("parent_trial_randomized", "transcriptomic_cohort_subjects", "baseline_available", "paired_week16", "paired_week52"),
  labels = c("Randomized", "Molecular\ncohort", "Baseline", "Paired\nW16", "Paired\nW52")
)]
flow_long[, treatment_label := factor(
  treatment,
  levels = c("Placebo", "Tabalumab_Q2W", "Tabalumab_Q4W"),
  labels = c("Placebo", "Q2W", "Q4W")
)]
flow_long[, randomized_fraction := participants / parent_trial_randomized]
flow_long[, trial := factor(trial, levels = c("ILLUMINATE-1", "ILLUMINATE-2"))]
p9a <- ggplot(flow_long, aes(stage, treatment_label, fill = randomized_fraction)) +
  geom_tile(colour = "white", linewidth = 0.9) +
  geom_text(aes(label = participants), size = 3.2, fontface = "bold") +
  facet_wrap(~trial, nrow = 1) +
  scale_fill_gradient(low = "#F3EFE8", high = "#287D78", limits = c(0, 1), labels = label_percent(accuracy = 1), name = "Fraction of\nrandomized arm") +
  labs(title = "Parent-trial randomization and longitudinal molecular availability", x = NULL, y = NULL) +
  theme_article(9.2) +
  theme(panel.grid = element_blank(), legend.position = "right")

component_heat <- copy(ancova_components)
component_heat[, component_label := factor(unname(component_labels[coordinate]), levels = unname(component_labels[coord_names]))]
component_heat[, contrast := paste(trial, unname(dose_labels[dose]), ifelse(time == "week16", "W16", "W52"), sep = " | ")]
component_heat[, contrast := factor(contrast, levels = rev(unique(contrast)))]
component_limit <- max(abs(component_heat$estimate), na.rm = TRUE)
p9b <- ggplot(component_heat, aes(component_label, contrast, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.75) +
  geom_tile(data = component_heat[fdr < 0.05], fill = NA, colour = "#151515", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%+.2f", estimate)), size = 2.9) +
  scale_fill_gradient2(low = "#2B6CB0", mid = "#F7F7F5", high = "#C53030", midpoint = 0,
                       limits = c(-component_limit, component_limit), name = "Drug-minus-\nplacebo effect") +
  labs(title = "Three-program treatment components", x = NULL, y = NULL) +
  theme_article(8.8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 24, hjust = 1), legend.position = "right")

moderation_trial <- read_tsv(moderation_trial_path)[, .(trial, time, dose, estimate, conf_low, conf_high, fdr)]
moderation_meta <- read_tsv(moderation_meta_path)[, .(trial = "Fixed-effect pooled", time, dose, estimate, conf_low, conf_high, fdr)]
moderation <- rbind(moderation_trial, moderation_meta)
moderation[, comparison := paste(unname(dose_labels[dose]), ifelse(time == "week16", "W16", "W52"), sep = " | ")]
moderation_levels <- rev(c("Q2W | W16", "Q4W | W16", "Q2W | W52", "Q4W | W52"))
moderation[, comparison := factor(comparison, levels = moderation_levels)]
moderation[, trial := factor(trial, levels = c("ILLUMINATE-1", "ILLUMINATE-2", "Fixed-effect pooled"))]
moderation_offsets <- c("ILLUMINATE-1" = 0.20, "ILLUMINATE-2" = 0, "Fixed-effect pooled" = -0.20)
moderation[, y_position := as.numeric(comparison) + unname(moderation_offsets[as.character(trial)])]
moderation_colors <- c(trial_colors, "Fixed-effect pooled" = "#252A31")
p9d <- ggplot(moderation, aes(colour = trial)) +
  geom_vline(xintercept = 0, colour = "#AEB4BA", linetype = "22", linewidth = 0.55) +
  geom_segment(aes(x = conf_low, xend = conf_high, y = y_position, yend = y_position), linewidth = 0.72) +
  geom_point(aes(x = estimate, y = y_position), size = 2.5) +
  scale_y_continuous(breaks = seq_along(moderation_levels), labels = moderation_levels, expand = expansion(add = c(0.42, 0.42))) +
  scale_colour_manual(values = moderation_colors, name = NULL) +
  labs(title = "Placebo-adjusted baseline B-coordinate interaction", subtitle = "Follow-up ANCOVA; effect per 1-SD higher baseline coordinate", x = "Treatment × baseline interaction", y = NULL) +
  theme_article(8.8) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

abundance_samples <- read_tsv(abundance_samples_path)[time == "baseline"]
abundance_samples[, trial := factor(trial, levels = c("ILLUMINATE-1", "ILLUMINATE-2"))]
abundance_correlations <- read_tsv(abundance_correlations_path)[
  subset == "baseline" & proxy == "quantiseq_b_fraction_z" & method == "spearman"
]
abundance_correlations[, trial := factor(trial, levels = c("ILLUMINATE-1", "ILLUMINATE-2"))]
abundance_correlations[, label := sprintf("Spearman ρ = %.3f\nn = %d", correlation, n)]
p9e <- ggplot(abundance_samples, aes(quantiseq_b_fraction_z, official_b_score)) +
  geom_point_rast(colour = "#2F7D67", alpha = 0.42, size = 1.25, raster.dpi = 450) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#1F5D4E",
              fill = "#9CC7BA", linewidth = 0.7, alpha = 0.22) +
  geom_label(data = abundance_correlations, aes(x = -Inf, y = Inf, label = label),
             inherit.aes = FALSE, hjust = -0.08, vjust = 1.15, size = 3.0,
             fill = "white", colour = "#202124", linewidth = 0.25) +
  facet_wrap(~trial, nrow = 1) +
  labs(title = "B-coordinate tracks estimated B-cell abundance",
       x = "quanTIseq B-cell fraction (standardized)", y = "B-cell-associated coordinate") +
  theme_article(9.0) +
  theme(legend.position = "none")

abundance_models <- read_tsv(abundance_models_path)[
  model %in% c("unadjusted_change", "quantiseq_adjusted_change", "marker30_adjusted_change")
]
abundance_models[, method_label := factor(
  model,
  levels = c("unadjusted_change", "quantiseq_adjusted_change", "marker30_adjusted_change"),
  labels = c("Unadjusted", "quanTIseq-adjusted", "Marker-proxy-adjusted")
)]
abundance_models[, contrast_label := paste(
  trial, unname(dose_labels[dose]), ifelse(time == "week16", "W16", "W52"), sep = " | "
)]
contrast_levels <- rev(c(
  "ILLUMINATE-1 | Q2W | W16", "ILLUMINATE-1 | Q4W | W16",
  "ILLUMINATE-1 | Q2W | W52", "ILLUMINATE-1 | Q4W | W52",
  "ILLUMINATE-2 | Q2W | W16", "ILLUMINATE-2 | Q4W | W16",
  "ILLUMINATE-2 | Q2W | W52", "ILLUMINATE-2 | Q4W | W52"
))
abundance_models[, contrast_label := factor(contrast_label, levels = contrast_levels)]
method_offsets <- c("Unadjusted" = 0.22, "quanTIseq-adjusted" = 0,
                    "Marker-proxy-adjusted" = -0.22)
abundance_models[, y_position := as.numeric(contrast_label) + unname(method_offsets[as.character(method_label)])]
abundance_models[, significance := factor(
  ifelse(fdr < 0.05, "FDR < 0.05", "FDR ≥ 0.05"), levels = c("FDR < 0.05", "FDR ≥ 0.05")
)]
model_colors <- c(
  "Unadjusted" = "#6B7280", "quanTIseq-adjusted" = "#2F7D67", "Marker-proxy-adjusted" = "#D18A27"
)
p9f <- ggplot(abundance_models, aes(colour = method_label)) +
  geom_vline(xintercept = 0, colour = "#AEB4BA", linetype = "22", linewidth = 0.55) +
  geom_segment(aes(x = bootstrap_conf_low, xend = bootstrap_conf_high,
                   y = y_position, yend = y_position), linewidth = 0.65) +
  geom_point(aes(x = estimate, y = y_position, shape = significance), size = 2.45, stroke = 0.85) +
  scale_y_continuous(breaks = seq_along(contrast_levels), labels = contrast_levels,
                     expand = expansion(add = c(0.45, 0.45))) +
  scale_colour_manual(values = model_colors, name = NULL) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR ≥ 0.05" = 1), name = NULL) +
  labs(title = "Treatment effects before and after abundance adjustment",
       x = "Drug-minus-placebo change in B-cell coordinate", y = NULL) +
  theme_article(8.5) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
        panel.grid.major.y = element_blank())

small <- theme(text=element_text(size=7),plot.title=element_text(size=7.5,face='bold',margin=margin(b=3)),
  plot.subtitle=element_blank(),axis.title=element_text(size=7),axis.text=element_text(size=6.5),
  legend.text=element_text(size=6.5),legend.title=element_text(size=6.5),
  strip.text=element_text(size=7),plot.margin=margin(3,3,3,3),
  legend.key.size=grid::unit(2.3,'mm'),legend.spacing.x=grid::unit(1,'mm'))
shrink_text <- function(p,size=2.2) {
  for(i in seq_along(p$layers)) if(inherits(p$layers[[i]]$geom,'GeomText') || inherits(p$layers[[i]]$geom,'GeomLabel')) p$layers[[i]]$aes_params$size<-size
  p+small
}
p9a <- shrink_text(p9a) + labs(title='Randomization and molecular sample availability',fill='Fraction') +
  scale_x_discrete(labels=c('Randomized','Molecular','Baseline','W16','W52')) +
  theme(axis.text.x=element_text(angle=30,hjust=1,size=6),legend.position='right') +
  guides(fill=guide_colourbar(barheight=grid::unit(15,'mm'),barwidth=grid::unit(2,'mm')))
p9b <- shrink_text(p9b,2.0) + labs(title='ANCOVA treatment effects',fill='Effect') +
  scale_y_discrete(labels=function(z) gsub('ILLUMINATE-','I',gsub(' \\| ',' ',z))) +
  scale_x_discrete(labels=c('Neutrophil','B cell','IFN')) +
  theme(axis.text.x=element_text(angle=0,hjust=.5,size=6),axis.text.y=element_text(size=6),legend.position='bottom') +
  guides(fill=guide_colourbar(barwidth=grid::unit(17,'mm'),barheight=grid::unit(2,'mm'),title.position='top'))
p9d <- p9d + small + labs(title='Baseline B-coordinate interaction',x='Treatment x baseline') +
  scale_colour_manual(values=moderation_colors,labels=c('I1','I2','Pooled'),name=NULL) +
  theme(legend.position='bottom',legend.text=element_text(size=6),axis.text.y=element_text(size=6))


p9e$layers <- p9e$layers[!vapply(p9e$layers,function(z) inherits(z$geom,'GeomSmooth'),logical(1))]
for(i in seq_along(p9e$layers)) if(inherits(p9e$layers[[i]]$geom,'GeomLabel')) {
  p9e$layers[[i]]$data[,label:=sprintf('Spearman rho = %.3f\nn = %d',correlation,n)]
}
p9e <- shrink_text(p9e,2.2) + labs(title='Baseline B-coordinate and B-cell abundance',
  x='quanTIseq B-cell fraction (standardized)',y='B-cell coordinate')
p9f <- p9f + small + labs(title='B-coordinate effects and abundance adjustment',
  x='Drug-minus-placebo change') +
  scale_y_continuous(breaks=seq_along(contrast_levels),labels=gsub('ILLUMINATE-','I',gsub(' \\| ',' ',contrast_levels)),
    expand=expansion(add=c(.45,.45))) +
  scale_colour_manual(values=model_colors,labels=c('Unadjusted','quanTIseq','Marker proxy'),name=NULL) +
  scale_shape_manual(values=c(16,1),labels=c('FDR < 0.05','FDR >= 0.05'),name=NULL) +
  theme(axis.text.y=element_text(size=6.5),legend.text=element_text(size=6.5),legend.box='horizontal')
fig <- p9a / (p9b | p9d) / p9e / p9f + plot_layout(heights=c(.85,1.15,.95,1.15))
fig <- fig + plot_annotation(tag_levels='A') & theme(plot.tag=element_text(family=font_family,size=9,face='bold'),plot.tag.position=c(0,.995))
ggsave(file.path(stage,'work/Figure8.pdf'),fig,width=5.12,height=6.5,units='in',device=cairo_pdf,family=font_family,bg='white')
ggsave(file.path(stage,'work/Figure8.png'),fig,width=5.12,height=6.5,units='in',device=ragg::agg_png,dpi=450,bg='white')
fwrite(data.table(panel=c('A','B','C','D','E'),source=c(flow_path,ancova_components_path,
 paste(moderation_trial_path,moderation_meta_path,sep=' | '),paste(abundance_samples_path,abundance_correlations_path,sep=' | '),abundance_models_path)),
 file.path(stage,'figures/Figure8_source_mapping.tsv'),sep='\t')
message('Figure8 rendered: A/B retained, original C removed, D/E/F become C/D/E.')
