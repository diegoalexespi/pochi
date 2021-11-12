#' Plot avereage feature expression per sample across groups in a Seurat object
#'
#' @param object Seurat object
#' @param features Which features to use for grouping cells
#' @param assay Location of features (typically a gene module assay, not RNA...)
#' @param split.by Which meta.data slot to use for splitting the groups
#' @param replicate.by Which meta.data slot to use as split.by replicates
#' @param draw_paths Whether or not to draw paths in the abundance plots between
#' replicates if they exist
#' @param perform.stat.test Whether or not to perform a statistical test between
#' the splits
#' @param test.choice Which choice of test to use (supports t.test of wilcox)
#' @param ncols Number of columns to plot the groups as
#' @param paired Whether the t.test should be paired if using t.test
#' @param set_minimum_to_zero Whether to set minimum abundance on plots to 0
#' globally
#' @param point_colors Vector of colors to color points as
#' @param p.adjust.method Which method to use for p-value adjustment for mulitple
#' testing
#' @param filter.p.above Plot p-values only below this filter
#' @param step.increase Tweaking the space between brackets for p-values vertically
#' @param title_size Text size of each group plot's title
#' @param p_val_size Text size for p-values
#' @param plot_type Either box or violin plot
#' @param sina_shift Whether to use ggforce::geom_sina to pretty-shift the points
#' @param x_lab The label on the x axis
#' @param y_lab The label on the y axis
#' @param point_size Size of the points plotted
#' @param same_y_limit Whether to keep the upper y limit constant across all
#' groups.
#' @details Plots the average expression of specific features in the Seurat
#' object across a split.by variable, using replicate.by as replicates for each
#' split.by condition. Typically used for gene module testing. Not recommended
#' for RNA or single-cell expression levels that are not Gaussian.
#'
#' @return Returns a patchwork plot of separate ggplots per features specified
#'
#' @importFrom rlang %||%
#'
#' @export
ModulePlot <- function(seurat_object,
                       assay = "AUC",
                       features = rownames(seurat_object[[assay]]),
                       split.by,
                       replicate.by,
                       draw_paths = FALSE,
                       perform.stat.test = TRUE,
                       test.choice = "t.test",
                       ncols = 3,
                       paired = FALSE,
                       set_minimum_to_zero = FALSE,
                       fill_colors = NULL,
                       point_colors = NULL,
                       p.adjust.method = "BH",
                       filter.p.above = 1,
                       step.increase = 0.15,
                       title_size = 10,
                       p_val_size = 2.5,
                       spacing_above = 1.3,
                       plot_type = "box",
                       sina_shift = TRUE,
                       x_lab = NULL,
                       y_lab = NULL,
                       same_y_limit = FALSE){

  x_lab <- if(is.null(x_lab)) split.by else x_lab
  y_lab <- if(is.null(y_lab)) "avg_expression" else y_lab

  my_splitter <- c(split.by, replicate.by)

  #nb the Seurat AverageExpression function returns colnames separated by a "_"
  avg_module_exp <- AverageExpression(seurat_object, assays = assay, group.by = my_splitter)[[1]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "module") %>%
    tidyr::pivot_longer(-module, values_to = "avg_score") %>%
    tidyr::separate(name, into = my_splitter, sep = "_") %>%
    tidyr::complete(nesting(!!sym(replicate.by), !!sym(split.by)), fill = list(avg_score = 0, n = 0)) %>%
    dplyr::mutate(split_by = !!sym(split.by), replicate_by = !!sym(replicate.by))

  facet.by.values <- features
  split.by.values <- unique(avg_module_exp$split_by)

  if(is.null(point_colors)){
    point_colors <- "black"
  }
  if(is.null(fill_colors)){
    fill_colors <- viridis::viridis_pal()(length(split.by.values))
  }

  my_plot_list <- lapply(seq_along(facet.by.values), function(i){
    target_group <- facet.by.values[i]
    target_group_scores <- avg_module_exp %>%
      dplyr::filter(module == target_group)
    if(perform.stat.test){
      t_test_results <- ggpubr::compare_means(avg_score~split_by,
                                              target_group_scores,
                                              paired = FALSE,
                                              p.adjust.method = p.adjust.method,
                                              method = test.choice)
    }
    #start ggplot
    g <- ggplot(target_group_scores, aes(x = split_by,
                                         y = avg_score,
                                         fill = split_by,
                                         group = split_by))+
      scale_fill_manual(values = fill_colors)+
      theme(title = element_text(size = title_size),
            panel.background = element_rect(fill = "white"),
            axis.line.x.bottom = element_line(color = 'black'),
            axis.line.y.left   = element_line(color = 'black'))+
      ggtitle(target_group)+
      xlab(x_lab)+
      NoLegend()

    #plot either violin or boxplot
    if(plot_type == "violin"){
      g <- g + geom_violin()
    } else if (plot_type == "box") {
      g <- g + geom_boxplot(outlier.shape = NA)
    } else {
      stop("plot_type must be one of violin or box")
    }


    #plot points as sina points or regular ggplot2 points
    if(sina_shift){
      g <- g + ggforce::geom_sina(color = point_colors)
    } else {
      g <- g + ggplot2::geom_point(color = point_colors)
    }

    #draw lines between replicates if applicable
    if(draw_paths){
      g <- g + geom_line(inherit.aes = FALSE, aes(x = split_by, y = avg_score, group = replicate_by))
    }

    #set plot minimum to 0 if needed (sometimes helpful)
    if(set_minimum_to_zero){
      ylim_min <- 0
    } else {
      ylim_min <- NA
    }

    #add brackets to statistical testing
    if(perform.stat.test){

      #get appropriate column and at the same time the prefix used
      p_column <- if(p.adjust.method == "none") "p" else "p.adj"

      #filtering if necessary
      t_test_results <- dplyr::filter(t_test_results, !!sym(p_column) <= filter.p.above)

      if(nrow(t_test_results) > 0){
        num_brackets <- nrow(t_test_results)
        max_score <- max(target_group_scores$avg_score)
        my_bracket_floor <- 1.05*max_score
        my_plot_ceiling <- max_score + my_bracket_floor*step.increase*num_brackets
        g <- g + ggpubr::geom_bracket(inherit.aes = FALSE,
                                      label.size = p_val_size,
                                      data = t_test_results,
                                      aes(xmin = group1,
                                          xmax = group2,
                                          label = paste0(p_column, ": ",signif(!!sym(p_column), 2))),
                                      y.position = my_bracket_floor,
                                      step.increase = step.increase)+
          scale_y_continuous(limits = c(ylim_min,my_plot_ceiling))
      } else {
        g <- g + scale_y_continuous(limits = c(ylim_min,NA))
      }

      #if no plotting of p values, do:
    } else {
      g <- g + scale_y_continuous(limits = c(ylim_min,NA))
    }

  })

  #get lowest y and highest y limit
  y_plot_limits <- lapply(my_plot_list, function(ggp) layer_scales(ggp)$y$get_limits())
  y_plot_limit_high <- max(sapply(y_plot_limits, "[[", 2))
  y_plot_limit_low <- min(sapply(y_plot_limits, "[[", 1))

  final_plot <- patchwork::wrap_plots(plotlist = my_plot_list, ncol = ncols)
  if(same_y_limit){
    final_plot <- final_plot &
      scale_y_continuous(limits = c(y_plot_limit_low, y_plot_limit_high), name = y_lab)
  } else {
    final_plot <- final_plot &
      scale_y_continuous(name = y_lab)
  }
  final_plot

}


