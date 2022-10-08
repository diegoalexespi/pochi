#' Plot cell abundances across groups in a Seurat object
#'
#' @param object Seurat object
#' @param group.by Which meta.data slot to use for grouping the cells
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
#' @param selected.groups Optional. If not NULL, will only plot the groups
#' specified in this argument.
#' @param rotated_axis Whether to apply Seurat::RotatedAxis to the plots
#' @details Plots the abundances of specific groups in the Seurat object
#' across a split.by variable, using replicate.by as replicates for each
#' split.by condition
#'
#' @return Returns a patchwork plot of separate ggplots per group.by variable
#'
#' @importFrom rlang %||%
#'
#' @export
AbundancePlot <- function(object,
                          group.by,
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
                          plot_type = "box",
                          sina_shift = TRUE,
                          x_lab = NULL,
                          y_lab = NULL,
                          point_size = 1,
                          same_y_limit = FALSE,
                          selected.groups = NULL,
                          rotated_axis = FALSE){

  x_lab <- if(is.null(x_lab)) split.by else x_lab
  y_lab <- if(is.null(y_lab)) "percentage" else y_lab

  seurat_metadata <- object@meta.data
  seurat_metadata_filtered <- seurat_metadata %>%
    dplyr::group_by(!!sym(replicate.by), !!sym(split.by), !!sym(group.by)) %>%
    dplyr::tally() %>%
    dplyr::mutate(Frequency = n/sum(n)) %>%
    dplyr::ungroup() %>%
    droplevels() %>%
    tidyr::complete(tidyr::nesting(!!sym(replicate.by), !!sym(split.by)), !!sym(group.by), fill = list(Frequency = 0, n = 0))

  split.by.values <- unique(seurat_metadata_filtered[[split.by]])

  frequencies <- lapply(seq_along(split.by.values), function(i){
    seurat_metadata_filtered %>%
      dplyr::filter(!!sym(split.by) == split.by.values[i]) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!sym(replicate.by), !!sym(group.by), Frequency) %>%
      dplyr::mutate(split_by = split.by.values[i], replicate_by = !!sym(replicate.by), group_by = !!sym(group.by))
  }) %>% do.call(rbind, .)

  if(!is.null(selected.groups)){
    frequencies <- frequencies %>% dplyr::filter(group_by %in% selected.groups)
    if(nrow(frequencies) < 1){
      stop("Must select at least 1 group if using selected.groups argument")
    }
  }

  group.by.names <- sort(unique(frequencies[["group_by"]]))

  if(is.null(point_colors)){
    point_colors <- "black"
  }
  if(is.null(fill_colors)){
    fill_colors <- viridis::viridis_pal()(length(split.by.values))
  }

  my_plot_list <- lapply(seq_along(group.by.names), function(i){
    target_group <- group.by.names[i]
    target_group_frequencies <- frequencies %>%
      dplyr::filter(group_by == target_group)
    if(perform.stat.test){
      t_test_results <- ggpubr::compare_means(Frequency~split_by,
                                              target_group_frequencies,
                                              paired = paired,
                                              p.adjust.method = p.adjust.method,
                                              method = test.choice)
    }

    #start ggplot
    g <- ggplot2::ggplot(target_group_frequencies, aes(x = split_by,
                                                       y = Frequency,
                                                       fill = split_by,
                                                       group = split_by))+
      ggplot2::scale_fill_manual(values = fill_colors)+
      ggplot2::theme(title = ggplot2::element_text(size = title_size),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     axis.line.x.bottom = ggplot2::element_line(color = 'black'),
                     axis.line.y.left   = ggplot2::element_line(color = 'black'))+
      ggplot2::ggtitle(target_group)+
      ggplot2::xlab(x_lab)+
      Seurat::NoLegend()

    #plot either violin or boxplot
    if(plot_type == "violin"){
      g <- g + ggplot2::geom_violin()
    } else if (plot_type == "box") {
      g <- g + ggplot2::geom_boxplot(outlier.shape = NA)
    } else {
      stop("plot_type must be one of violin or box")
    }


    #plot points as sina points or regular ggplot2 points
    if(sina_shift){
      g <- g + ggforce::geom_sina(color = point_colors, size = point_size)
    } else {
      g <- g + ggplot2::geom_point(color = point_colors, size = point_size)
    }

    #draw lines between replicates if applicable
    if(draw_paths){
      g <- g + ggplot2::geom_line(inherit.aes = FALSE, aes(x = split_by, y = Frequency, group = replicate_by))
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
        max_frequency <- max(target_group_frequencies$Frequency)
        my_bracket_floor <- 1.05*max_frequency
        my_plot_ceiling <- max_frequency + my_bracket_floor*step.increase*num_brackets
        g <- g + ggpubr::geom_bracket(inherit.aes = FALSE,
                                      label.size = p_val_size,
                                      data = t_test_results,
                                      aes(xmin = group1,
                                          xmax = group2,
                                          label = paste0(p_column, ": ",signif(!!sym(p_column), 2))),
                                      y.position = my_bracket_floor,
                                      step.increase = step.increase)+
          ggplot2::scale_y_continuous(limits = c(ylim_min,my_plot_ceiling), labels = function(x) paste0(x*100, "%"), name = y_lab)
      } else {
        g <- g + ggplot2::scale_y_continuous(limits = c(ylim_min,NA), labels = function(x) paste0(x*100, "%"), name = y_lab)
      }

      #if no plotting of p values, do:
    } else {
      g <- g + ggplot2::scale_y_continuous(limits = c(ylim_min,NA), labels = function(x) paste0(x*100, "%"), name = y_lab)
    }

  })

  #get lowest y and highest y limit
  y_plot_limits <- lapply(my_plot_list, function(ggp) layer_scales(ggp)$y$get_limits())
  y_plot_limit_high <- max(sapply(y_plot_limits, "[[", 2))
  y_plot_limit_low <- min(sapply(y_plot_limits, "[[", 1))


  final_plot <- patchwork::wrap_plots(plotlist = my_plot_list, ncol = ncols)
  if(same_y_limit){
    final_plot <- final_plot &
      ggplot2::scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits = c(y_plot_limit_low, y_plot_limit_high), name = y_lab)
  }
  if(rotated_axis){
    final_plot <- final_plot &
      Seurat::RotatedAxis()
  }
  final_plot
}
