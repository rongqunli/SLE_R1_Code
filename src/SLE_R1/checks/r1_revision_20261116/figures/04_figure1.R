if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")


options(stringsAsFactors=FALSE,scipen=99)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales)})
stage <- 'F:/program/Hei/SLE_R1/checks/r1_revision_20261116'
dir.create(file.path(stage,'work'),recursive=TRUE,showWarnings=FALSE)
base <- 'F:/program/Hei/SLE_v1'
font_family <- 'Times New Roman'
target_modules <- c('magenta','greenyellow','cyan')
theme_article <- function(base_size=7) theme_classic(base_size=base_size,base_family=font_family) +
  theme(plot.title=element_text(face='bold',size=7.5,margin=margin(b=3)),
    axis.text=element_text(colour='#202124',size=6.5),axis.title=element_text(size=7),
    panel.grid.major.y=element_line(colour='#EEEEEE',linewidth=.2),
    strip.background=element_rect(fill='#F3F4F5',colour='#D5D8DC',linewidth=.3),
    strip.text=element_text(size=6.5,face='bold'),legend.position='bottom',
    legend.text=element_text(size=6),legend.title=element_text(size=6.5),plot.margin=margin(3,3,3,3))
read_result <- function(name) fread(file.path(base,'results/3.1/06_wgcna',name))
soft_diagnostics <- read_result('soft_threshold_diagnostics.tsv')
associations <- read_result('module_trait_associations.tsv')
stability <- read_result('module_stability.tsv')
annotation_terms <- fread(file.path(base,'Figures/tables/selected_annotation_terms_for_Figure1.tsv'))
soft_diagnostics[, power_label := as.character(power)]
p1b_left <- ggplot(soft_diagnostics, aes(power, SFT.R.sq)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "#777777", linewidth = 0.55) +
  geom_line(colour = "#B33B35", linewidth = 0.55) +
  geom_text(aes(label = power_label, fontface = ifelse(selected, "bold", "plain")),
            colour = "#B33B35", size = 3.2, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(2, 12, 2)) +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "Soft-threshold power", y = "Signed scale-free fit (R²)") +
  theme_article(9.5)
p1b_right <- ggplot(soft_diagnostics, aes(power, mean.k.)) +
  geom_line(colour = "#B33B35", linewidth = 0.55) +
  geom_text(aes(label = power_label, fontface = ifelse(selected, "bold", "plain")),
            colour = "#B33B35", size = 3.2, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(2, 12, 2)) +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  labs(x = "Soft-threshold power", y = "Mean connectivity") +
  theme_article(9.5)
p1b <- (p1b_left | p1b_right) +
  plot_annotation(title = "Soft-threshold selection") &
  theme(plot.title = element_text(family = font_family, face = "bold", size = 12.5,
                                  margin = margin(b = 6, l = 14)))
module_traits <- associations[model_family %chin% c(
  "activity_composition_adjusted", "renal_composition_adjusted",
  "c3_composition_adjusted", "c4_composition_adjusted"
)]
module_traits[, aligned_t := fifelse(grepl("^c[34]_", model_family), -t_statistic, t_statistic)]
module_traits[, axis_x := fifelse(model_family == "activity_composition_adjusted", 1,
                          fifelse(model_family == "renal_composition_adjusted", 2,
                            fifelse(model_family == "c3_composition_adjusted", 3, 4)))]
axis_labels <- c("SLEDAI", "Renal activity", "Low C3", "Low C4")
other_order <- module_traits[!module %chin% target_modules,
                             .(selection_score = max(aligned_t)), by = module][order(-selection_score), module]
top_to_bottom <- c(target_modules, other_order)
module_traits[, module_display := factor(module, levels = rev(top_to_bottom))]
module_traits[, selected_module := module %chin% target_modules]
module_traits[, cell_label := paste0(
  "t=", sprintf("%+.1f", aligned_t), "\nq=",
  ifelse(fdr < 0.001, formatC(fdr, format = "e", digits = 1), sprintf("%.3f", fdr))
)]
strip_data <- unique(module_traits[, .(module, module_display)])
p1d <- ggplot(module_traits, aes(axis_x, module_display, fill = pmax(-8, pmin(8, aligned_t)))) +
  geom_tile(width = 0.94, height = 0.92, colour = "white", linewidth = 0.5) +
  geom_tile(data = module_traits[selected_module == TRUE],
            width = 0.94, height = 0.92, fill = NA, colour = "#111111", linewidth = 0.8) +
  geom_text(aes(label = cell_label), size = 3.0, lineheight = 0.87) +
  geom_tile(data = strip_data, aes(x = 0.42, y = module_display),
            width = 0.16, height = 0.92, fill = strip_data$module,
            colour = "#444444", linewidth = 0.2, inherit.aes = FALSE) +
  scale_fill_gradient2(low = "#3566A8", mid = "white", high = "#C84A46",
                       midpoint = 0, limits = c(-8, 8), name = "Activity-aligned t") +
  scale_x_continuous(breaks = 1:4, labels = axis_labels,
                     limits = c(0.25, 4.55), expand = expansion(mult = 0)) +
  scale_y_discrete(labels = function(x) paste0("ME", x), expand = expansion(add = 0.35)) +
  labs(title = "Adjusted module-trait landscape", x = NULL, y = NULL) +
  theme_article(9.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold"),
        legend.position = "right")
stability_plot <- stability[module %chin% target_modules]
stability_plot[, module_display := factor(paste0("ME", module),
                                           levels = paste0("ME", rev(target_modules)))]
p1e <- ggplot(stability_plot,
              aes(median_kme_rank_correlation, module_display, colour = module)) +
  geom_segment(aes(x = q025_kme_rank_correlation,
                   xend = q975_kme_rank_correlation, yend = module_display),
               linewidth = 1.2) +
  geom_point(size = 3.2) +
  geom_text(aes(x = q975_kme_rank_correlation + 0.0007, label = paste0("n=", n_genes)),
            hjust = 0, colour = "#222222", size = 3.2) +
  scale_colour_manual(values = setNames(target_modules, target_modules), guide = "none") +
  scale_x_continuous(limits = c(0.968, 0.997), breaks = seq(0.97, 0.995, 0.005)) +
  labs(title = "Module membership stability",
       x = "Patient-bootstrap kME rank correlation", y = NULL) +
  theme_article(10)
small <- theme_article(7)
shrink_text <- function(p,size=2.0) {
  for(i in seq_along(p$layers)) if(inherits(p$layers[[i]]$geom,'GeomText')) p$layers[[i]]$aes_params$size<-size
  p+small
}
p1b_left <- shrink_text(p1b_left) + labs(title='Scale-free fit',y=expression(R^2),x='Power') +
  scale_y_continuous(breaks=c(0,.4,.8))
p1b_right <- shrink_text(p1b_right) + labs(title='Mean connectivity',y='Mean connectivity',x='Power')
p1b <- p1b_left | p1b_right


p1d$data$cell_label <- sprintf('%+.1f',p1d$data$aligned_t)
p1d <- shrink_text(p1d,size=1.9) +
  theme(legend.position='right',axis.text.y=element_text(size=6.3),axis.text.x=element_text(size=6.5)) +
  guides(fill=guide_colourbar(title='Aligned t',barheight=grid::unit(18,'mm'),barwidth=grid::unit(2,'mm')))
p1e <- shrink_text(p1e,size=2.0) + labs(x='Patient-bootstrap kME rank correlation') +
  scale_x_continuous(limits=c(.968,.998),breaks=c(.97,.98,.99))

annotation_terms[,module_display:=factor(module,levels=rev(target_modules))]
annotation_terms[,term_id:=factor(paste(module,pathway,sep='___'),levels=rev(paste(module,pathway,sep='___')))]
term_labels <- setNames(annotation_terms$display_term,as.character(annotation_terms$term_id))
p1f <- ggplot(annotation_terms,aes(fold_enrichment,term_id,size=overlap,fill=neglog10_fdr,colour=module)) +
  geom_point(shape=21,stroke=.5,alpha=.96) +
  geom_text(aes(x=108,label=fdr_label),hjust=1,size=1.8,colour='#202124') +
  facet_grid(module_display~.,scales='free_y',space='free_y') +
  scale_y_discrete(labels=function(z) unname(term_labels[z])) +
  scale_colour_manual(values=setNames(target_modules,target_modules),guide='none') +
  scale_fill_gradient(low='#E8F0F6',high='#7A1F4D',name='-log10(FDR)') +
  scale_size_continuous(range=c(1.2,3),breaks=c(20,40),name='Overlap') +
  scale_x_continuous(limits=c(0,110),breaks=c(0,20,40,60,80),expand=expansion(mult=c(.01,0))) +
  guides(fill=guide_colourbar(barwidth=grid::unit(14,'mm'),barheight=grid::unit(2,'mm')),
    size=guide_legend(override.aes=list(fill='white'))) +
  labs(title='GO and KEGG annotation of selected modules',x='Fold enrichment',y=NULL) + small +
  theme(axis.text.y=element_text(size=5.7),strip.text.y=element_text(angle=0,size=5.7),
    legend.key.size=grid::unit(2,'mm'),legend.box='horizontal')

fig <- wrap_elements(full=p1b) / wrap_elements(full=p1d) /
  wrap_elements(full=p1e) / wrap_elements(full=p1f) +
  plot_layout(heights=c(1.0,1.9,.9,2.7))
fig <- fig + plot_annotation(tag_levels='A') &
  theme(plot.tag=element_text(family=font_family,size=9,face='bold'),plot.tag.position=c(0,.995))
ggsave(file.path(stage,'work/Figure1.pdf'),fig,width=5.12,height=6.5,units='in',device=cairo_pdf,family=font_family,bg='white')
ggsave(file.path(stage,'work/Figure1.png'),fig,width=5.12,height=6.5,units='in',device=ragg::agg_png,dpi=450,bg='white')
fwrite(annotation_terms,file.path(stage,'figures/Figure1_enrichment_plotted_values.tsv'),sep='\t')
message('Figure1: four retained panels, saved results, dendrogram removed.')
