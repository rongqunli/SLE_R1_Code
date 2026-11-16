if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")


options(stringsAsFactors = FALSE, scipen = 99)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

project_root <- "F:/program/Hei/SLE_v1"
output_root <- file.path(project_root, "Figures", "Figure3")
record_root <- file.path(
  project_root, "results", "reviewer_response", "02_clinical_landscape_validation"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

font_family <- "Times New Roman"
state_labels <- c(
  HEALTHY_CONTROL = "Healthy control",
  S01_NEUTROPHIL_HIGH = "Neutrophil-high",
  S02_IMMUNE_LOW = "Immune-low",
  S03_B_CELL_HIGH = "B-cell-high",
  S04_IFN_HIGH = "IFN-high"
)
state_colors <- c(
  "Healthy control" = "#3A9D78", "Neutrophil-high" = "#D18A27",
  "Immune-low" = "#707A86", "B-cell-high" = "#2F7D67", "IFN-high" = "#A23B3B"
)
model_labels <- c(
  binary_high_vs_low = "Binary high vs low",
  distance_from_immune_low = "Distance from immune-low",
  four_group = "Four groups",
  continuous_three_program = "Three continuous programs"
)
model_colors <- c(
  "Binary high vs low" = "#8C8C8C", "Distance from immune-low" = "#D18A27",
  "Four groups" = "#4C78A8", "Three continuous programs" = "#2F7D67"
)
trait_labels <- c(
  SLEDAI = "SLEDAI", C3 = "C3", C4 = "C4", anti_dsDNA = "anti-dsDNA"
)

theme_article <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, colour = "#202124"),
      plot.title = element_text(face = "bold", size = rel(1.16), margin = margin(b = 7)),
      plot.subtitle = element_blank(), plot.caption = element_blank(),
      axis.text = element_text(colour = "#202124"),
      strip.background = element_rect(fill = "#F3F4F5", colour = "#D5D8DC", linewidth = 0.45),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.3),
      panel.grid.minor = element_blank(), legend.position = "bottom",
      legend.title = element_blank(), plot.margin = margin(9, 11, 9, 9)
    )
}
theme_set(theme_article())
update_geom_defaults("text", list(family = font_family, colour = "#202124"))

save_plot <- function(plot, stem, width, height) {
  ggsave(
    file.path(output_root, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, family = font_family, bg = "white"
  )
  ggsave(
    file.path(output_root, paste0(stem, ".png")), plot,
    width = width, height = height, units = "in", dpi = 320,
    device = ragg::agg_png, background = "white"
  )
}

pca_path <- file.path(
  project_root, "results", "exploratory_program_landscape", "02_continuous_landscape",
  "adult_program_landscape_pca.tsv"
)
pca_data <- fread(pca_path)
pca_data[, state_label := factor(unname(state_labels[state]), levels = unname(state_labels))]
p3a <- ggplot(pca_data, aes(PC1, PC2, colour = state_label)) +
  geom_point(alpha = 0.38, size = 1.05) +
  stat_ellipse(
    data = pca_data[state != "HEALTHY_CONTROL"], aes(group = state_label),
    type = "norm", level = 0.80, linewidth = 0.75, show.legend = FALSE
  ) +
  scale_colour_manual(values = state_colors, name = NULL) +
  labs(title = "Continuous three-program molecular landscape", x = "Landscape PC1", y = "Landscape PC2") +
  theme_article(10)
save_plot(p3a, "Figure3A_continuous_program_landscape", 8.2, 6.2)

profiles_path <- file.path(
  project_root, "results", "exploratory_program_landscape", "02_continuous_landscape",
  "adult_state_program_profiles.tsv"
)
profiles_wide <- fread(profiles_path)
profiles <- melt(
  profiles_wide,
  id.vars = c("state", "n", "mean_landscape_distance", "median_landscape_distance"),
  measure.vars = grep("^coord_", names(profiles_wide), value = TRUE),
  variable.name = "program", value.name = "mean_coordinate"
)
profiles[, state_label := factor(unname(state_labels[state]), levels = rev(unname(state_labels)))]
profiles[, program_label := factor(c(
  coord_CYAN_NEUTROPHIL_IMMATURE = "Immature-neutrophil",
  coord_GREENYELLOW_B_CELL_COMPLEMENT = "B-cell/BCR",
  coord_MAGENTA_IFN_ANTIVIRAL = "IFN/antiviral"
)[program], levels = c("IFN/antiviral", "B-cell/BCR", "Immature-neutrophil"))]
p3b <- ggplot(profiles, aes(program_label, state_label, fill = mean_coordinate)) +
  geom_tile(colour = "white", linewidth = 0.9) +
  geom_text(aes(label = sprintf("%+.2f", mean_coordinate)), size = 3.35) +
  scale_fill_gradient2(
    low = "#3566A8", mid = "white", high = "#C84A46", midpoint = 0,
    name = "Mean coordinate"
  ) +
  labs(title = "Frozen nearest-centroid region profiles", x = NULL, y = NULL) +
  theme_article(9.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")
save_plot(p3b, "Figure3B_reference_state_program_profiles", 7.4, 4.8)

bootstrap_path <- file.path(record_root, "bootstrap_optimism_corrected_metrics.tsv")
cv_path <- file.path(record_root, "repeated_cross_validation_summary.tsv")
bootstrap <- fread(bootstrap_path)[model %chin% names(model_labels)]
cross_validation <- fread(cv_path)[model %chin% names(model_labels)]
internal_bootstrap <- bootstrap[, .(
  trait, model,
  validation_method = "Bootstrap optimism-corrected",
  estimate = optimism_corrected_incremental_r_squared,
  conf_low = incremental_r2_bootstrap_q025,
  conf_high = incremental_r2_bootstrap_q975
)]
internal_cv <- cross_validation[, .(
  trait, model,
  validation_method = "Repeated 10-fold cross-validation",
  estimate = repeated_cv_incremental_r2_mean,
  conf_low = repeated_cv_incremental_r2_q025,
  conf_high = repeated_cv_incremental_r2_q975
)]
internal <- rbindlist(list(internal_bootstrap, internal_cv))
internal[, `:=`(
  trait_label = factor(unname(trait_labels[trait]), levels = unname(trait_labels)),
  model_label = factor(unname(model_labels[model]), levels = unname(model_labels)),
  validation_method = factor(
    validation_method,
    levels = c("Bootstrap optimism-corrected", "Repeated 10-fold cross-validation")
  )
)]
p3c <- ggplot(internal, aes(trait_label, estimate, colour = model_label, group = model_label)) +
  geom_hline(yintercept = 0, colour = "#A8ABB0", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = conf_low, ymax = conf_high), width = 0.10,
    position = position_dodge(width = 0.66), linewidth = 0.55
  ) +
  geom_point(position = position_dodge(width = 0.66), size = 2.4) +
  facet_wrap(~validation_method, nrow = 1) +
  scale_colour_manual(values = model_colors, name = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.08, 0.12))) +
  labs(title = "Internal validation of clinical variance captured", x = NULL, y = "Incremental R-squared") +
  theme_article(9.3) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
save_plot(p3c, "Figure3C_internal_clinical_validation", 13.5, 5.8)

heldout_path <- file.path(record_root, "trial_linked_heldout_validation_metrics.tsv")
heldout <- fread(heldout_path)[model %chin% c("four_group", "continuous_three_program")]
validation_labels <- c(
  ILLUMINATE1_validation = "ILLUMINATE-1 linked",
  ILLUMINATE2_validation = "ILLUMINATE-2 linked",
  pooled_trial_linked_validation = "Pooled trial linked"
)
heldout[, `:=`(
  trait_label = factor(unname(trait_labels[trait]), levels = unname(trait_labels)),
  model_label = factor(unname(model_labels[model]), levels = unname(model_labels[c("four_group", "continuous_three_program")])),
  validation_label = factor(unname(validation_labels[validation_set]), levels = unname(validation_labels))
)]
p3d <- ggplot(heldout, aes(trait_label, heldout_incremental_r_squared, colour = model_label, group = model_label)) +
  geom_hline(yintercept = 0, colour = "#A8ABB0", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = heldout_incremental_r2_q025, ymax = heldout_incremental_r2_q975),
    width = 0.10, position = position_dodge(width = 0.46), linewidth = 0.55
  ) +
  geom_point(position = position_dodge(width = 0.46), size = 2.5) +
  facet_wrap(~validation_label, nrow = 1) +
  scale_colour_manual(values = model_colors, name = NULL, drop = TRUE) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.08, 0.12))) +
  labs(title = "Fixed-model performance in trial-linked held-out participants", x = NULL, y = "Held-out incremental R-squared") +
  theme_article(9.3) +
  theme(axis.text.x = element_text(angle = 24, hjust = 1), legend.position = "bottom")
save_plot(p3d, "Figure3D_trial_linked_heldout_validation", 14.5, 5.5)

top_row <- (p3a | p3b) + plot_layout(widths = c(1.08, 0.92))
figure3 <- top_row / p3c / p3d + plot_layout(heights = c(1.02, 0.88, 0.88))
figure3_tagged <- figure3 + plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(family = font_family, face = "bold", size = 15),
    plot.tag.position = c(0.002, 0.998)
  )
save_plot(figure3_tagged, "Figure3_combined", 16.0, 17.0)

manifest <- data.table(
  panel = c("A", "B", "C", "D", "combined"),
  source = c(
    pca_path, profiles_path, paste(bootstrap_path, cv_path, sep = " | "),
    heldout_path, paste(pca_path, profiles_path, bootstrap_path, cv_path, heldout_path, sep = " | ")
  )
)
fwrite(manifest, file.path(record_root, "Figure3_validation_manifest.tsv"), sep = "\t")
message("Figure 3 clinical validation rebuild completed: ", output_root)

stage <- 'F:/program/Hei/SLE_R1/checks/r1_revision_20261116/work'
dir.create(stage,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(stage,'Figure3.pdf'),figure3_tagged,width=5.12,height=5.8,units='in',device=cairo_pdf,family=font_family,bg='white')
