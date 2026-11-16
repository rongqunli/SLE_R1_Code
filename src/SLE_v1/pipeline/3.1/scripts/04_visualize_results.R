if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(c("data.table", "ggplot2", "patchwork"),
                 "final figures and Chinese result summary")

script_path <- script_file()
output_dir <- file.path(cfg$paths$results, "04_visualization")
summary_file <- file.path(output_dir, "GSE65391_program_validation结果总结.md")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  file.path(cfg$paths$results, "00_input_audit", "input_audit_summary.tsv"),
  file.path(cfg$paths$results, "01_expression_qc", "pca_coordinates.tsv"),
  file.path(cfg$paths$results, "01_expression_qc", "sample_qc_metrics.tsv"),
  file.path(cfg$paths$results, "02_program_scores", "program_gene_matching.tsv"),
  file.path(cfg$paths$results, "03_longitudinal", "longitudinal_analysis_data.tsv"),
  file.path(cfg$paths$results, "03_longitudinal", "longitudinal_model_results.tsv"),
  file.path(cfg$paths$results, "03_longitudinal", "primary_leave_one_patient_out.tsv"),
  file.path(cfg$paths$results, "03_longitudinal", "anti_dsDNA_format_audit.tsv")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Missing visualization input(s): ", paste(missing_inputs, collapse = "; "))
}

figure_names <- c(
  "fig01_study_flow.png",
  "fig02_pca_qc.png",
  "fig03_within_patient_sledai.png",
  "fig04_model_forest.png",
  "fig05_leave_one_patient_out.png",
  "fig06_program_coverage_composition.png",
  "fig07_secondary_associations.png"
)
pdf_name <- "GSE65391_visualization_book.pdf"
manifest_file <- file.path(output_dir, "figure_manifest.tsv")
outputs <- c(file.path(output_dir, c(figure_names, pdf_name)), manifest_file, summary_file)
signature <- build_step_signature(
  input_files = c(inputs, cfg$files$programs),
  code_files = c(script_path, file.path(PIPELINE_ROOT, "scripts", "_common.R"),
                 file.path(PIPELINE_ROOT, "config", "config.R")),
  extra = list(version = "seven_figure_book_v1")
)
if (!has_flag("--force") && step_is_complete(output_dir, signature, outputs)) {
  log_message("Step 04 visualization is already complete.")
  quit(save = "no", status = 0L)
}

read_tsv <- function(path) data.table::fread(path, sep = "\t", na.strings = "NA")
audit <- read_tsv(inputs[[1]])
pca <- read_tsv(inputs[[2]])
qc <- read_tsv(inputs[[3]])
matching <- read_tsv(inputs[[4]])
analysis <- read_tsv(inputs[[5]])
models <- read_tsv(inputs[[6]])
loo <- read_tsv(inputs[[7]])
dsdna <- read_tsv(inputs[[8]])
program_config <- read_tsv(cfg$files$programs)

format_number <- function(x) formatC(as.numeric(x), format = "g", digits = 8)
format_p <- function(x) {
  value <- as.numeric(x)
  ifelse(value < 0.001, formatC(value, format = "e", digits = 3),
         formatC(value, format = "f", digits = 4))
}
theme_pipeline <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(color = "grey20")
    )
}
save_png <- function(plot, filename, width, height) {
  ggplot2::ggsave(
    file.path(output_dir, filename), plot = plot, width = width, height = height,
    units = "in", dpi = 220, bg = "white"
  )
}
metric_value <- function(name) {
  value <- audit[metric == name, value]
  if (length(value) != 1L) stop("Missing audit metric: ", name)
  as.numeric(value)
}


flow <- data.table::data.table(
  stage = factor(
    c("GEO matrix", "Included arrays", "SLE visits", "SLE patients"),
    levels = c("GEO matrix", "Included arrays", "SLE visits", "SLE patients")
  ),
  x = 1:4,
  count = c(
    metric_value("matrix_arrays"), metric_value("included_arrays"),
    metric_value("included_sle_visits"), metric_value("sle_patients")
  ),
  detail = c(
    "all normalized arrays", "24 technical replicates excluded",
    "48 healthy arrays reserved for QC", "patient-level inference"
  ),
  fill = c("#4477AA", "#44AA99", "#CC6677", "#AA3377")
)
p1 <- ggplot2::ggplot(flow, ggplot2::aes(x, 1)) +
  ggplot2::geom_segment(
    data = flow[x < 4], ggplot2::aes(x = x + 0.36, xend = x + 0.64, y = 1, yend = 1),
    arrow = grid::arrow(length = grid::unit(0.18, "inches")), linewidth = 0.8,
    color = "grey45", inherit.aes = FALSE
  ) +
  ggplot2::geom_tile(ggplot2::aes(fill = fill), width = 0.7, height = 0.62,
                     color = "white", linewidth = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(stage, "\n", count)),
                     color = "white", fontface = "bold", size = 4.3, lineheight = 1.15) +
  ggplot2::geom_text(ggplot2::aes(y = 0.57, label = detail), size = 3.1,
                     color = "grey25", lineheight = 1.05) +
  ggplot2::scale_fill_identity() +
  ggplot2::coord_cartesian(xlim = c(0.55, 4.45), ylim = c(0.35, 1.4), clip = "off") +
  ggplot2::labs(
    title = "GSE65391 study flow",
    subtitle = "Technical arrays are audited separately from biological visits and patients"
  ) +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(color = "grey30"),
    plot.margin = ggplot2::margin(15, 20, 30, 20)
  )
save_png(p1, figure_names[[1]], 13, 4.8)


pca_qc <- merge(
  pca,
  qc[, .(sample_id, median_correlation, robust_pc_distance, outlier_flag)],
  by = "sample_id", all.x = TRUE, sort = FALSE
)
pca_qc[, outlier_label := factor(
  ifelse(outlier_flag, "QC flagged", "Retained"),
  levels = c("Retained", "QC flagged")
)]
pca_qc[, disease_state := factor(disease_state, levels = c("Healthy", "SLE"))]
p2a <- ggplot2::ggplot(
  pca_qc,
  ggplot2::aes(PC1, PC2, color = disease_state, shape = outlier_label)
) +
  ggplot2::geom_point(alpha = 0.72, size = 2.1) +
  ggplot2::scale_color_manual(values = c("Healthy" = "#4477AA", "SLE" = "#CC6677")) +
  ggplot2::scale_shape_manual(values = c("Retained" = 16, "QC flagged" = 4)) +
  ggplot2::labs(title = "Expression PCA", x = "PC1", y = "PC2",
                color = "Disease state", shape = "QC status") +
  theme_pipeline()
p2b <- ggplot2::ggplot(
  pca_qc,
  ggplot2::aes(median_correlation, robust_pc_distance, color = outlier_label)
) +
  ggplot2::geom_point(alpha = 0.65, size = 2) +
  ggplot2::scale_color_manual(values = c("Retained" = "grey62", "QC flagged" = "#D55E00")) +
  ggplot2::labs(title = "Pre-defined outlier metrics", x = "Median sample correlation",
                y = "Robust PCA distance", color = "QC status") +
  theme_pipeline()
p2 <- (p2a | p2b) + patchwork::plot_annotation(
  title = "Expression-level quality control",
  subtitle = "Six flagged SLE arrays are retained in the primary analysis and excluded only in sensitivity analysis"
)
save_png(p2, figure_names[[2]], 14, 6.5)


analysis[, cytotoxicity_within := CYTOTOXICITY_SHARED - mean(CYTOTOXICITY_SHARED),
         by = patient_id]
analysis[, qc_label := factor(
  ifelse(qc_outlier_flag, "QC flagged", "Primary set"),
  levels = c("Primary set", "QC flagged")
)]
p3 <- ggplot2::ggplot(
  analysis,
  ggplot2::aes(SLEDAI_within, cytotoxicity_within, color = qc_label)
) +
  ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.4) +
  ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.4) +
  ggplot2::geom_point(alpha = 0.42, size = 1.7) +
  ggplot2::geom_smooth(
    data = analysis[qc_outlier_flag == FALSE], method = "lm", formula = y ~ x,
    se = TRUE, color = "#0072B2", fill = "#56B4E9", linewidth = 1
  ) +
  ggplot2::scale_color_manual(values = c("Primary set" = "grey55", "QC flagged" = "#D55E00")) +
  ggplot2::labs(
    title = "Within-patient disease activity and shared cytotoxicity",
    subtitle = "Points are visits; both axes remove each patient's mean for descriptive display",
    x = "Within-patient SLEDAI deviation",
    y = "Within-patient CYTOTOXICITY_SHARED deviation", color = NULL
  ) + theme_pipeline(11)
save_png(p3, figure_names[[3]], 10.5, 7)


family_labels <- c(
  sle_activity_base = "Primary",
  sle_activity_repeated_patients = "Repeated patients",
  sle_activity_qc_sensitivity = "QC sensitivity",
  sle_activity_composition_adjusted = "Composition adjusted",
  sle_activity_clinical_adjusted = "Clinical adjusted"
)
forest <- models[
  program %in% c("CYTOTOXICITY_SHARED", "CD8_GZMK_STATE") &
    model_family %in% names(family_labels) & term == "SLEDAI_within"
]
forest[, family_label := factor(
  family_labels[model_family], levels = rev(unname(family_labels))
)]
forest[, estimate_scaled := estimate * 1000]
forest[, low_scaled := conf_low * 1000]
forest[, high_scaled := conf_high * 1000]
forest[, supported := factor(
  ifelse(conf_low > 0 | conf_high < 0, "95% CI excludes zero", "95% CI crosses zero"),
  levels = c("95% CI excludes zero", "95% CI crosses zero")
)]
p4 <- ggplot2::ggplot(
  forest,
  ggplot2::aes(estimate_scaled, family_label, color = supported)
) +
  ggplot2::geom_vline(xintercept = 0, color = "grey45", linetype = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = low_scaled, xmax = high_scaled),
    orientation = "y", width = 0.18, linewidth = 0.75
  ) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_wrap(~program, ncol = 1, scales = "free_x") +
  ggplot2::scale_color_manual(values = c(
    "95% CI excludes zero" = "#D55E00", "95% CI crosses zero" = "#0072B2"
  )) +
  ggplot2::labs(
    title = "Frozen-program longitudinal models",
    subtitle = "Effect per one-point within-patient SLEDAI increase; estimates multiplied by 1,000",
    x = "Model estimate x 1,000", y = NULL, color = NULL
  ) + theme_pipeline(10)
save_png(p4, figure_names[[4]], 11, 8)


data.table::setorder(loo, estimate)
loo[, exclusion_rank := seq_len(.N)]
primary <- models[
  model_family == "sle_activity_base" & program == cfg$longitudinal$primary_program &
    term == "SLEDAI_within"
]
if (nrow(primary) != 1L) stop("Expected one primary model row.")
p5 <- ggplot2::ggplot(loo, ggplot2::aes(exclusion_rank, estimate * 1000)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey45", linetype = 2) +
  ggplot2::geom_hline(yintercept = primary$estimate * 1000, color = "#D55E00", linewidth = 0.9) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = conf_low * 1000, ymax = conf_high * 1000),
    fill = "#56B4E9", alpha = 0.18
  ) +
  ggplot2::geom_point(color = "#0072B2", size = 1.3, alpha = 0.75) +
  ggplot2::labs(
    title = "Leave-one-patient-out stability",
    subtitle = "All 158 re-fits retain a negative primary estimate; orange line is the full-cohort estimate",
    x = "Excluded patient, ranked by estimate", y = "Estimate x 1,000"
  ) + theme_pipeline(11)
save_png(p5, figure_names[[5]], 11, 6.5)


coverage <- matching[, .(
  matched_genes = sum(matched), configured_genes = .N
), by = .(program, category)]
coverage[, coverage_fraction := matched_genes / configured_genes]
coverage[, program := factor(program, levels = rev(program_config$program))]
p6a <- ggplot2::ggplot(coverage, ggplot2::aes(coverage_fraction, program, fill = category)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(matched_genes, "/", configured_genes)),
    hjust = 1.12, color = "white", fontface = "bold", size = 3.2
  ) +
  ggplot2::scale_fill_manual(values = c("biological" = "#AA3377", "composition" = "#44AA99")) +
  ggplot2::scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1)) +
  ggplot2::labs(title = "Frozen program coverage", x = "Matched fraction", y = NULL,
                fill = "Program class") + theme_pipeline(10)
composition_compare <- models[
  program %in% c("CYTOTOXICITY_SHARED", "CD8_GZMK_STATE") &
    model_family %in% c("sle_activity_base", "sle_activity_composition_adjusted") &
    term == "SLEDAI_within"
]
composition_compare[, model_label := factor(
  ifelse(model_family == "sle_activity_base", "Primary", "Composition adjusted"),
  levels = c("Composition adjusted", "Primary")
)]
p6b <- ggplot2::ggplot(
  composition_compare,
  ggplot2::aes(estimate * 1000, model_label, color = program)
) +
  ggplot2::geom_vline(xintercept = 0, color = "grey45", linetype = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = conf_low * 1000, xmax = conf_high * 1000),
    orientation = "y", width = 0.16, linewidth = 0.75
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = c(
    "CYTOTOXICITY_SHARED" = "#D55E00", "CD8_GZMK_STATE" = "#0072B2"
  )) +
  ggplot2::labs(
    title = "Composition changes the primary interpretation",
    x = "SLEDAI effect x 1,000", y = NULL, color = "Frozen program"
  ) + theme_pipeline(10)
p6 <- (p6a | p6b) + patchwork::plot_annotation(
  title = "Program coverage and composition sensitivity"
)
save_png(p6, figure_names[[6]], 14, 7)


secondary_labels <- c(
  renal_activity_association = "Renal activity",
  c3_association = "C3 within-patient",
  c4_association = "C4 within-patient"
)
secondary <- models[
  program %in% c("CYTOTOXICITY_SHARED", "CD8_GZMK_STATE") &
    model_family %in% names(secondary_labels)
]
secondary[, outcome_label := factor(
  secondary_labels[model_family], levels = rev(unname(secondary_labels))
)]
secondary[, supported := factor(
  ifelse(fdr < cfg$longitudinal$fdr, "FDR < 0.05", "FDR >= 0.05"),
  levels = c("FDR < 0.05", "FDR >= 0.05")
)]
p7 <- ggplot2::ggplot(secondary, ggplot2::aes(estimate, outcome_label, color = supported)) +
  ggplot2::geom_vline(xintercept = 0, color = "grey45", linetype = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = conf_low, xmax = conf_high),
    orientation = "y", width = 0.17, linewidth = 0.75
  ) +
  ggplot2::geom_point(size = 2.9) +
  ggplot2::facet_wrap(~program, ncol = 1, scales = "free_x") +
  ggplot2::scale_color_manual(values = c("FDR < 0.05" = "#009E73", "FDR >= 0.05" = "grey55")) +
  ggplot2::labs(
    title = "Secondary renal and complement associations",
    subtitle = "C3/C4 source units are unavailable; effect sizes are not comparable across outcomes",
    x = "Model estimate", y = NULL, color = NULL
  ) + theme_pipeline(10)
save_png(p7, figure_names[[7]], 11, 7.5)

plots <- list(p1, p2, p3, p4, p5, p6, p7)
pdf_file <- file.path(output_dir, pdf_name)
grDevices::cairo_pdf(pdf_file, width = 14, height = 9, onefile = TRUE)
for (plot_object in plots) print(plot_object)
grDevices::dev.off()


primary <- models[
  model_family == "sle_activity_base" & program == "CYTOTOXICITY_SHARED" &
    term == "SLEDAI_within"
]
repeated <- models[
  model_family == "sle_activity_repeated_patients" & program == "CYTOTOXICITY_SHARED"
]
qc_sensitive <- models[
  model_family == "sle_activity_qc_sensitivity" & program == "CYTOTOXICITY_SHARED"
]
composition <- models[
  model_family == "sle_activity_composition_adjusted" & program == "CYTOTOXICITY_SHARED"
]
clinical <- models[
  model_family == "sle_activity_clinical_adjusted" & program == "CYTOTOXICITY_SHARED"
]
cd8_base <- models[
  model_family == "sle_activity_base" & program == "CD8_GZMK_STATE"
]
cd8_composition <- models[
  model_family == "sle_activity_composition_adjusted" & program == "CD8_GZMK_STATE"
]
selected_rows <- list(primary, repeated, qc_sensitive, composition, clinical, cd8_base, cd8_composition)
if (any(vapply(selected_rows, nrow, integer(1)) != 1L)) stop("Missing key summary model row.")

model_table_line <- function(label, row) {
  paste0(
    "| ", label, " | `", format_number(row$estimate), "` | [`",
    format_number(row$conf_low), "`, `", format_number(row$conf_high), "`] | `",
    format_p(row$p_value), "` | `", format_p(row$fdr), "` | ",
    row$n_patients, " | ", row$n_visits, " |"
  )
}
coverage_shared <- coverage[program == "CYTOTOXICITY_SHARED"]
coverage_cd8 <- coverage[program == "CD8_GZMK_STATE"]
secondary_bio <- secondary[order(model_family, program)]
secondary_lines <- vapply(seq_len(nrow(secondary_bio)), function(i) {
  row <- secondary_bio[i]
  paste0(
    "| ", secondary_labels[[row$model_family]], " | ", row$program, " | `",
    format_number(row$estimate), "` | `", format_p(row$p_value), "` | `",
    format_p(row$fdr), "` |"
  )
}, character(1))
loo_negative <- sum(loo$estimate < 0)
summary_lines <- c(
  "# GSE65391_pipeline结果总结",
  "",
  paste0("> 总结日期：", format(Sys.Date(), "%Y-%m-%d"), "  "),
  "> 状态：步骤00–04已完成；结论按患者内纵向效应解释",
  "",
  "## 1. 阶段性结论",
  "",
  paste0(
    "GSE65391未复现GSE254176中细胞毒/CD8 GZMK程序随疾病活动升高的正向信号。",
    "在924次SLE访视、158名患者的主模型中，`CYTOTOXICITY_SHARED`的患者内SLEDAI效应为`",
    format_number(primary$estimate), "`，95% CI [`", format_number(primary$conf_low), "`, `",
    format_number(primary$conf_high), "`]，P=`", format_number(primary$p_value), "`，FDR=`",
    format_number(primary$fdr), "`；方向为负，与GSE254176方向相反。"
  ),
  "",
  paste0(
    "该负向效应在仅重复访视患者、排除6个QC标记访视、临床调整以及",
    loo_negative, "/", nrow(loo), "次留一患者分析中保持负向；但在细胞组成调整后变为`",
    format_number(composition$estimate), "`（95% CI [`", format_number(composition$conf_low), "`, `",
    format_number(composition$conf_high), "`]，P=`", format_number(composition$p_value),
    "`），接近零且不显著。"
  ),
  "",
  paste0(
    "因此最稳妥的判断是：本全血队列不支持一个可独立于细胞组成的正向CD8 GZMK/细胞毒疾病活动轴。",
    "未调整模型中的负相关很可能依赖NK/T细胞等细胞组成，不能解释为CD8细胞内在状态或因果机制。"
  ),
  "",
  "## 2. 数据与分析对象",
  "",
  "| 项目 | 数量 | 分析角色 |",
  "|---|---:|---|",
  paste0("| GEO标准化数组 | ", metric_value("matrix_arrays"), " | 输入审计 |"),
  paste0("| 纳入数组 | ", metric_value("included_arrays"), " | 表达质控与程序评分 |"),
  paste0("| 排除技术重复 | ", metric_value("technical_replicates"), " | 不作为独立生物学观察 |"),
  paste0("| SLE访视 | ", metric_value("included_sle_visits"), " | 纵向模型 |"),
  paste0("| SLE患者 | ", metric_value("sle_patients"), " | 主要推断单位 |"),
  paste0("| 健康数组 | ", metric_value("healthy_arrays"), " | 仅用于表达QC |"),
  "",
  "患者而不是数组或探针是统计单位。主要模型将SLEDAI拆分为患者内偏差和患者间均值，并使用患者随机截距。",
  "",
  "## 3. 主要与稳健性结果",
  "",
  "| 模型 | 估计值 | 95% CI | P值 | FDR | 患者 | 访视 |",
  "|---|---:|---:|---:|---:|---:|---:|",
  model_table_line("共享细胞毒：主模型", primary),
  model_table_line("共享细胞毒：仅重复访视患者", repeated),
  model_table_line("共享细胞毒：排除QC标记访视", qc_sensitive),
  model_table_line("共享细胞毒：组成调整", composition),
  model_table_line("共享细胞毒：临床调整", clinical),
  model_table_line("CD8 GZMK：主模型", cd8_base),
  model_table_line("CD8 GZMK：组成调整", cd8_composition),
  "",
  paste0(
    "留一患者分析共", nrow(loo), "次，", loo_negative, "次估计为负；估计范围为`",
    format_number(min(loo$estimate)), "`至`", format_number(max(loo$estimate)), "`。",
    "这排除了单一患者驱动主模型方向的可能，但不能消除细胞组成混杂。"
  ),
  "",
  "## 4. 冻结程序覆盖与次要结果",
  "",
  paste0(
    "主要`CYTOTOXICITY_SHARED`程序匹配", coverage_shared$matched_genes, "/",
    coverage_shared$configured_genes, "个基因；`CD8_GZMK_STATE`匹配",
    coverage_cd8$matched_genes, "/", coverage_cd8$configured_genes,
    "个基因。程序在查看GSE65391结局前冻结，未用本队列结果挑选基因或改变方向。"
  ),
  "",
  "| 次要结局 | 程序 | 估计值 | P值 | FDR |",
  "|---|---|---:|---:|---:|",
  secondary_lines,
  "",
  "肾脏活动、C3和C4属于预设次要分析。C3/C4来源单位缺失，因此只能解释关联方向和队列内相对变化，不能跨研究比较效应量。anti-dsDNA同时含数值、滴度、阴性文本和缺失值，本版仅完成格式审计，不合并为连续变量。",
  "",
  "## 5. 质量控制与解释边界",
  "",
  "- 6个预定义表达QC离群访视在主分析中保留，并在敏感性分析中排除；排除后患者内SLEDAI重新中心化。",
  "- GSE65391是全血芯片，程序分数混合细胞比例与细胞内状态，不能单独证明CD8 GZMK细胞数量、状态或TCR克隆发生变化。",
  "- 混合模型报告关联而不是因果；当前结果不能直接触发虚拟敲除、药物重定位或结构筛选。",
  "- 健康样本用于QC，不替代患者内疾病活动比较。",
  "",
  "## 6. 下一步决策",
  "",
  "预设门槛要求主要效应为正、置信区间支持且组成调整与留一患者方向稳定。本队列没有达到该门槛，因此不应基于GSE65391把该轴升级到GSE121239或下游靶点筛选。应冻结本次未复现结论，并回到跨队列设计层面区分全血组成变化与细胞类型特异状态。",
  "",
  "## 7. 图表与机器结果",
  "",
  "全部独立PNG、合并PDF和图表清单位于`results/04_visualization/`；模型明细位于`results/03_longitudinal/longitudinal_model_results.tsv`。"
)
temporary_summary <- paste0(summary_file, ".tmp-", Sys.getpid())
on.exit(unlink(temporary_summary, force = TRUE), add = TRUE)
writeLines(summary_lines, temporary_summary, useBytes = TRUE)
if (file.exists(summary_file)) unlink(summary_file, force = TRUE)
if (!file.rename(temporary_summary, summary_file)) stop("Could not write summary file.")

figure_paths <- file.path(output_dir, c(figure_names, pdf_name))
manifest <- data.table::data.table(
  figure = c(sprintf("Figure %d", 1:7), "Combined book"),
  file = basename(figure_paths),
  type = c(rep("png", 7), "pdf"),
  title = c(
    "Study flow", "Expression PCA and QC", "Within-patient SLEDAI association",
    "Frozen-program model forest", "Leave-one-patient-out stability",
    "Program coverage and composition sensitivity", "Secondary associations",
    "Seven-page visualization book"
  ),
  interpretation_status = c(
    "design", "quality_control", "descriptive", "primary_and_sensitivity",
    "robustness", "coverage_and_robustness", "secondary", "combined"
  ),
  bytes = as.numeric(file.info(figure_paths)$size)
)
safe_write_tsv(manifest, manifest_file)
write_session_info(file.path(cfg$paths$logs, "04_visualization_sessionInfo.txt"))
mark_step_complete(output_dir, signature, outputs)
log_message("Visualization and summary complete: 7 PNG figures, one PDF book, and Chinese summary.")
