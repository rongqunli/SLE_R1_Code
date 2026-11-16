render_figure5_atlas_panels <- function(metadata, palette, output_root, record_root, python, python_helper) {
  for (package in c('data.table', 'ggplot2', 'ggrastr', 'ggrepel', 'ggunchull', 'circlize', 'MASS', 'png', 'ragg')) {
    if (!requireNamespace(package, quietly = TRUE)) stop('Missing package: ', package)
  }
  atlas <- data.table::copy(data.table::as.data.table(metadata))
  palette <- data.table::copy(data.table::as.data.table(palette))
  data.table::setorder(palette, order)
  needed <- c('cell_id', 'UMAP_1', 'UMAP_2', 'subject_id', 'condition', 'cell_type_analysis')
  stopifnot(all(needed %in% names(atlas)), !anyDuplicated(atlas$cell_id), all(is.finite(atlas$UMAP_1)), all(is.finite(atlas$UMAP_2)), setequal(unique(atlas$cell_type_analysis), palette$cell_type_analysis), setequal(unique(as.character(atlas$condition)), c('HD', 'SLE')))
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  records <- file.path(record_root, 'AB')
  dir.create(records, recursive = TRUE, showWarnings = FALSE)
  invoke <- function(arguments) {
    status <- system2(python, args = shQuote(c(python_helper, arguments)))
    if (!identical(status, 0L)) stop('Figure 5 PDF helper failed with status ', status)
    invisible(status)
  }
  atlas[, cell_type := factor(cell_type_analysis, levels = palette$cell_type_analysis)]
  colours <- setNames(palette$colour, palette$cell_type_analysis)
  a_labels <- aggregate(as.data.frame(atlas[, .(UMAP_1, UMAP_2)]), list(cell_type = atlas$cell_type), median)
  a_labels$label <- palette$short_label[match(a_labels$cell_type, palette$cell_type_analysis)]
  coordinate_min <- floor(min(atlas$UMAP_1, atlas$UMAP_2))
  lower <- coordinate_min - .1 * abs(coordinate_min)
  label_relative <- .1 * abs(lower)
  end <- lower + abs(.25 * lower)
  middle <- lower + abs(.25 * lower) / 2
  axes <- data.frame(x = c(lower, lower), xend = c(end, lower), y = c(lower, lower), yend = c(lower, end))
  axis_labels <- data.frame(x = c(middle, lower - label_relative), y = c(lower - label_relative, middle), label = c('UMAP1', 'UMAP2'), angle = c(0, 90))
  set.seed(20260718L)
  p <- ggplot2::ggplot(atlas, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggrastr::geom_point_rast(ggplot2::aes(colour = cell_type), size = .19, alpha = .68, raster.dpi = 450) +
    ggunchull::stat_unchull(ggplot2::aes(fill = cell_type, group = cell_type), delta = .2, alpha = .2, linewidth = 1.5, colour = 'grey50', linetype = 'dashed', show.legend = FALSE) +
    ggplot2::geom_segment(data = axes, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, linewidth = .45, colour = 'black', arrow = grid::arrow(length = grid::unit(.10, 'inches'), type = 'closed')) +
    ggplot2::geom_text(data = axis_labels, ggplot2::aes(x = x, y = y, label = label, angle = angle), inherit.aes = FALSE, family = 'Arial', fontface = 'italic', size = 3.4) +
    ggrepel::geom_text_repel(data = a_labels, ggplot2::aes(UMAP_1, UMAP_2, label = label), inherit.aes = FALSE, family = 'Arial', colour = 'black', size = 4.3, seed = 20260718L, box.padding = .50, point.padding = .18, force = 2, min.segment.length = 0, segment.colour = 'grey45', segment.size = .35, max.overlaps = Inf, max.iter = 20000, max.time = 5) +
    ggplot2::scale_colour_manual(values = colours, breaks = palette$cell_type_analysis, labels = palette$label, drop = FALSE, name = NULL) +
    ggplot2::scale_fill_manual(values = colours, drop = FALSE) +
    ggplot2::coord_equal(clip = 'off') +
    ggplot2::labs(title = 'PBMC atlas', x = NULL, y = NULL, tag = 'A') +
    ggplot2::theme_bw(base_family = 'Arial', base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.border = ggplot2::element_rect(fill = NA, colour = 'black', linewidth = .65), legend.position = 'right', legend.text = ggplot2::element_text(size = 10), legend.key.height = grid::unit(5, 'mm'), legend.key.width = grid::unit(5, 'mm'), plot.title = ggplot2::element_text(size = 14, face = 'bold', hjust = .5), plot.tag = ggplot2::element_text(size = 17, face = 'bold'), plot.tag.position = c(0, 1), plot.margin = ggplot2::margin(12, 14, 12, 12)) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1), ncol = 1))
  built <- ggplot2::ggplot_build(p)
  stopifnot(nrow(built$data[[1]]) == nrow(atlas), identical(as.numeric(built$data[[1]]$x), atlas$UMAP_1), identical(as.numeric(built$data[[1]]$y), atlas$UMAP_2), length(unique(built$data[[2]]$group)) == nrow(palette))
  ragg::agg_capture(width = 9, height = 7, units = 'in', res = 450)
  grob <- ggplot2::ggplot_gtable(built)
  grDevices::dev.off()
  a_pdf <- file.path(output_root, 'Figure5A_PBMC_atlas.pdf')
  ggplot2::ggsave(a_pdf, grob, width = 9, height = 7, units = 'in', device = grDevices::cairo_pdf, family = 'Arial', bg = 'white')
  ggplot2::ggsave(sub('pdf$', 'png', a_pdf), grob, width = 9, height = 7, units = 'in', device = ragg::agg_png, dpi = 450, bg = 'white')
  counts <- atlas[, .(cells = .N, donors = data.table::uniqueN(subject_id)), by = condition]
  data.table::setorder(counts, condition)
  center_x <- mean(range(atlas$UMAP_1))
  center_y <- mean(range(atlas$UMAP_2))
  radius <- sqrt((atlas$UMAP_1 - center_x)^2 + (atlas$UMAP_2 - center_y)^2)
  coordinate_scale <- .78 / max(radius)
  atlas[, `:=`(x = (UMAP_1 - center_x) * coordinate_scale, y = (UMAP_2 - center_y) * coordinate_scale)]
  stopifnot(max(abs(atlas$x / coordinate_scale + center_x - atlas$UMAP_1)) < 1e-12, max(abs(atlas$y / coordinate_scale + center_y - atlas$UMAP_2)) < 1e-12)
  atlas[, cell_type_analysis := factor(cell_type_analysis, levels = palette$cell_type_analysis)]
  atlas[, colour := palette$colour[match(cell_type_analysis, palette$cell_type_analysis)]]
  bandwidth <- c(MASS::bandwidth.nrd(atlas$x), MASS::bandwidth.nrd(atlas$y))
  contour_levels <- c(.2, .3)
  density <- lapply(c('HD', 'SLE'), function(group) {
    data <- atlas[condition == group]
    MASS::kde2d(data$x, data$y, h = bandwidth, n = 200, lims = c(-.88, .88, -.88, .88))
  })
  names(density) <- c('HD', 'SLE')
  stopifnot(identical(density$HD$x, density$SLE$x), identical(density$HD$y, density$SLE$y))
  centers <- atlas[, .(anchor_x = median(x), anchor_y = median(y), cells = .N), by = cell_type_analysis]
  centers[, label := palette$short_label[match(cell_type_analysis, palette$cell_type_analysis)]]
  data.table::setorder(centers, -cells)
  layout_labels <- function() {
    placed <- list()
    for (i in seq_len(nrow(centers))) {
      row <- centers[i]
      half_width <- strwidth(row$label, cex = .62) / 2 + .010
      half_height <- strheight('Mg', cex = .62) / 2 + .007
      radii <- c(0, .025, .05, .075, .10, .125, .15, .175)
      angles <- seq(0, 2 * pi, length.out = 25)[-25]
      candidates <- rbind(data.table::data.table(x = row$anchor_x, y = row$anchor_y), data.table::rbindlist(lapply(radii[-1], function(distance) data.table::data.table(x = row$anchor_x + distance * cos(angles), y = row$anchor_y + distance * sin(angles)))))
      selected <- NULL
      for (j in seq_len(nrow(candidates))) {
        candidate <- candidates[j]
        box <- c(candidate$x - half_width, candidate$x + half_width, candidate$y - half_height, candidate$y + half_height)
        inside <- max(sqrt(c(box[1], box[2])^2 + max(abs(box[3:4]))^2)) < .91
        overlap <- length(placed) > 0 && any(vapply(placed, function(old) box[1] < old$xmax && box[2] > old$xmin && box[3] < old$ymax && box[4] > old$ymin, logical(1)))
        if (inside && !overlap) {
          selected <- cbind(row, data.table::data.table(x = candidate$x, y = candidate$y, xmin = box[1], xmax = box[2], ymin = box[3], ymax = box[4]))
          break
        }
      }
      if (is.null(selected)) stop('Unable to place label: ', row$label)
      placed[[length(placed) + 1L]] <- selected
    }
    data.table::rbindlist(placed)
  }
  labels <- NULL
  layer_metadata <- list()
  layer_images <- list()
  draw_points_layer <- function(data, group) {
    usr <- par('usr')
    physical <- par('pin')
    if (is.null(layer_images[[group]])) {
      path <- file.path(records, paste0('points_', group, '_600dpi.png'))
      vector_path <- file.path(records, paste0('points_', group, '_cairo.pdf'))
      point_size <- par('ps')
      base_cex <- par('cex')
      previous_device <- grDevices::dev.cur()
      grDevices::cairo_pdf(vector_path, width = physical[1], height = physical[2], family = 'Arial', bg = 'transparent', pointsize = point_size)
      par(mar = rep(0, 4), xaxs = 'i', yaxs = 'i', bg = NA, cex = base_cex, ps = point_size)
      plot.new()
      plot.window(xlim = usr[1:2], ylim = usr[3:4], xaxs = 'i', yaxs = 'i')
      points(data$x, data$y, pch = 16, col = grDevices::adjustcolor(data$colour, alpha.f = .20), cex = .10)
      grDevices::dev.off()
      grDevices::dev.set(previous_device)
      invoke(c('render', vector_path, path, '--dpi', '600'))
      layer_images[[group]] <<- png::readPNG(path)
      layer_metadata[[group]] <<- data.table::data.table(group = group, cells = nrow(data), donors = data.table::uniqueN(data$subject_id), dpi = 600, width_in = physical[1], height_in = physical[2], width_pixels = dim(layer_images[[group]])[2], height_pixels = dim(layer_images[[group]])[1], xmin = usr[1], xmax = usr[2], ymin = usr[3], ymax = usr[4], point_cex = .1, point_alpha = .2)
    }
    saved <- layer_metadata[[group]]
    stopifnot(max(abs(usr - unlist(saved[, .(xmin, xmax, ymin, ymax)], use.names = FALSE))) < 1e-12, max(abs(physical - unlist(saved[, .(width_in, height_in)], use.names = FALSE))) < 1e-12)
    rasterImage(layer_images[[group]], saved$xmin, saved$ymin, saved$xmax, saved$ymax, interpolate = FALSE)
  }
  plot_group <- function(group) {
    data <- atlas[condition == group]
    circlize::circos.clear()
    circlize::circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(.01, 0), track.height = .012, gap.degree = rep(2, nrow(palette)), start.degree = 90, clock.wise = TRUE, canvas.xlim = c(-1.17, 1.17), canvas.ylim = c(-1.17, 1.17), points.overflow.warning = FALSE)
    circlize::circos.initialize(factors = factor(palette$cell_type_analysis, levels = palette$cell_type_analysis), xlim = cbind(rep(0, nrow(palette)), rep(1, nrow(palette))))
    circlize::circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = .012, panel.fun = function(x, y) {
      type <- circlize::get.cell.meta.data('sector.index')
      index <- match(type, palette$cell_type_analysis)
      circlize::circos.lines(c(0, 1), c(.5, .5), col = palette$colour[index], lwd = 3)
      circlize::circos.text(.5, 1 + circlize::mm_y(3.4), palette$short_label[index], cex = .68, col = '#202124', facing = 'bending.inside', niceFacing = TRUE)
    })
    draw_points_layer(data, group)
    contour(density[[group]], drawlabels = FALSE, levels = contour_levels, col = '#ae9c76', lwd = .65, add = TRUE)
    if (is.null(labels)) labels <<- layout_labels()
    moved <- sqrt((labels$x - labels$anchor_x)^2 + (labels$y - labels$anchor_y)^2) > .035
    segments(labels$anchor_x[moved], labels$anchor_y[moved], labels$x[moved], labels$y[moved], col = '#666666', lwd = .4)
    rect(labels$xmin + .003, labels$ymin + .003, labels$xmax - .003, labels$ymax - .003, col = grDevices::adjustcolor('white', alpha.f = .88), border = NA)
    text(labels$x, labels$y, labels$label, cex = .62, col = '#202124')
    current <- counts[condition == group]
    mtext(group, side = 3, line = 1.5, cex = .95, font = 2)
    mtext(sprintf('%d donors | %s cells', current$donors, format(current$cells, big.mark = ',', scientific = FALSE, trim = TRUE)), side = 3, line = .35, cex = .72)
    circlize::circos.clear()
  }
  b_pdf <- file.path(output_root, 'Figure5B_Circular_HD_SLE_for_AI.pdf')
  grDevices::cairo_pdf(b_pdf, width = 11.6, height = 6, family = 'Arial', bg = 'white', onefile = TRUE)
  par(family = 'Arial', mfrow = c(1, 2), mar = c(.4, .4, 3.1, .4), oma = c(.15, .1, .1, .1), xpd = NA, cex = 1)
  for (group in c('HD', 'SLE')) plot_group(group)
  grDevices::dev.off()
  invoke(c('render', b_pdf, sub('pdf$', 'png', b_pdf), '--dpi', '450', '--opaque'))
  for (group in c('HD', 'SLE')) {
    grDevices::cairo_pdf(file.path(records, paste0('Figure5B_Circular_', group, '_uncropped.pdf')), width = 11.6, height = 6, family = 'Arial', bg = 'transparent', onefile = TRUE)
    par(family = 'Arial', mfrow = c(1, 2), mar = c(.4, .4, 3.1, .4), oma = c(.15, .1, .1, .1), xpd = NA, cex = 1, bg = NA)
    plot_group(group)
    grDevices::dev.off()
  }
  invoke(c('crop-pair', file.path(records, 'Figure5B_Circular_HD_uncropped.pdf'), file.path(records, 'Figure5B_Circular_SLE_uncropped.pdf'), output_root))
  data.table::fwrite(counts, file.path(records, 'group_counts.tsv'), sep = '\t')
  data.table::fwrite(palette, file.path(records, 'celltype_label_colour_mapping.tsv'), sep = '\t')
  data.table::fwrite(labels, file.path(records, 'fixed_label_positions.tsv'), sep = '\t')
  data.table::fwrite(data.table::rbindlist(layer_metadata), file.path(records, 'point_layer_parameters.tsv'), sep = '\t')
  data.table::fwrite(data.table::data.table(parameter = c('center_x', 'center_y', 'shared_scale', 'max_transformed_radius', 'inverse_max_error_x', 'inverse_max_error_y', 'bandwidth_x', 'bandwidth_y', 'grid_n', 'grid_min', 'grid_max', 'contour_level_1', 'contour_level_2', 'point_cex', 'point_alpha', 'sectors', 'sector_width', 'gap_degrees'), value = c(center_x, center_y, coordinate_scale, max(sqrt(atlas$x^2 + atlas$y^2)), max(abs(atlas$x / coordinate_scale + center_x - atlas$UMAP_1)), max(abs(atlas$y / coordinate_scale + center_y - atlas$UMAP_2)), bandwidth, 200, -.88, .88, contour_levels, .1, .2, nrow(palette), 1, 2)), file.path(records, 'shared_plot_parameters.tsv'), sep = '\t')
  saveRDS(density, file.path(records, 'common_grid_density.rds'), compress = 'gzip')
  list(A = a_pdf, B = b_pdf)
}
